
//-ML:  Threshold for alpha1
scalar alpha1Lim = EL_alphaThreshold;   
//-ML:  Minimum number of cells to consider as a structure
label minCells = EL_cellThreshold;      

//-ML:  Reset volumeID 
volScalarField vofID_
(
    IOobject
    (
        "vofID",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("vofID", dimless, scalar(0.0)),
    zeroGradientFvPatchScalarField::typeName
);

//-ML:  Local-to-global cell reference
globalIndex globalNumbering (mesh.nCells());

//-ML:  Create the volumeID field
forAll(mesh.cells(), cellI)
{
    //-ML:  Convert local processor cellIDs to global cellIDs
    label globalCellID = globalNumbering.toGlobal(cellI);

    //-ML:  Cell is considered vapor and not yet marked
    if (alpha1[cellI] < alpha1Lim && mag(vofID_[cellI]) < 0.1)
    {

        label currentID = globalCellID+1;
        label startCell = cellI;
        //-ML:  Set volume ID
        vofID_[startCell] = currentID;

        //-ML:  Create dynamic list to store all neighbours
        DynamicList<label> cells(0);
        //-ML:  Get neighbours of startCell
        cells.append(mesh.cellCells()[startCell]);

        //-ML:  Start outer loop
        while (cells.size() > 0)
        {
            //-ML:  Create new dynamic list to store all new neighbours
            DynamicList<label> neighbours(0);

            //-ML:  Loop over all current neighbours
            forAll(cells, i)
            {
                label cellJ = cells[i];

                //-ML:  Check if neighbour cell is liquid and has no ID yet
                if (alpha1[cellJ] < alpha1Lim && mag(vofID_[cellJ] - currentID) > 0.1)
                {
                    //-ML:  Set volume ID
                    vofID_[cellJ] = currentID;

                    //-ML:  add neighbours to neighboursNeighbours list
                    neighbours.append(mesh.cellCells()[cellJ]);
                }
            }
            cells.clear();
            cells.transfer(neighbours);
        }
    }
}
vofID_.correctBoundaryConditions();


//-ML:  Correct for parallel processing
//-ML:  Number of corrections
label corr = Foam::log(scalar(Pstream::nProcs())) / Foam::log(2.0) + 1;
for(int i=0; i<=corr; i++)
{
    forAll(mesh.boundary(), patchI)
    {
        //-ML:  Current patch
        const fvPatch& pPatch = mesh.boundary()[patchI];

        //-ML:  check if patch is coupled, e.g. processor patch
        if (pPatch.coupled())
        {
            //-ML:  loop over all faces of processor patch
            forAll(pPatch, faceI)
            {
                //-ML:  get cell values
                const scalar ID_own = vofID_.boundaryField()[patchI].patchInternalField()()[faceI];
                const scalar ID_nei = vofID_.boundaryField()[patchI].patchNeighbourField()()[faceI];

                //-ML:  check, if volumeIDs on both sides are fluid
                if ( (ID_nei > 0.1) && (ID_own > 0.1) )
                {
                    //-ML:  check, if volumeIDs on both sides are different
                    if ( mag(ID_nei - ID_own) > 0.1 )
                    {
                        scalar minID = Foam::min(ID_own, ID_nei);
                        scalar maxID = Foam::max(ID_own, ID_nei);

                        scalar diff(maxID - minID);
                        volScalarField filtered = diff
                                                * pos( vofID_ - maxID + 0.1)
                                                * pos(-vofID_ + maxID + 0.1);
                        vofID_ -= filtered;
                    }
                }
            }
        }
    }
    vofID_.correctBoundaryConditions();
}

//-ML: Count number of continua and renumber vofID
{
    volScalarField tmp(vofID_);
    vofID_ = 0.0;
    scalar maxID = gMax(tmp);
    int iter = 1;
    while (maxID > 0.1)
    {
        volScalarField curVolID = pos(tmp - maxID + 0.1);
        vofID_ += iter*curVolID;
        tmp -= maxID*curVolID;
        maxID = gMax(tmp);
        iter++;
    }
}

//-ML: Calculate droplet volume, velocity, and position
//-ML: Create lists to store data
label maxID = floor(gMax(vofID_) + 0.5);
labelList noCells(maxID+1, label(0));
scalarList cellVolume(maxID+1, scalar(0));
scalarList volume(maxID+1, scalar(0));
vectorList position(maxID+1, vector(0,0,0));
vectorList velocity(maxID+1, vector(0,0,0));

//-ML: Loop over all cells and store corresponding data
forAll(mesh.cells(), cellI)
{
    if (vofID_[cellI] > 0.1)
    {
        //-ML: Get volume ID
        label volID = floor(vofID_[cellI] + 0.5);

        //-ML: Store data in corresponding list
        scalar oneMinusalpha1 = (1 - alpha1[cellI]);
        noCells[volID] += 1;
        cellVolume[volID] += mesh.V()[cellI];
        volume[volID] += oneMinusalpha1*mesh.V()[cellI];
        velocity[volID] += oneMinusalpha1*mesh.V()[cellI]*U[cellI];
        position[volID] += oneMinusalpha1*mesh.V()[cellI]*mesh.C()[cellI];
    }
}

//-ML: Account for parallel processing
reduce( noCells, sumOp<labelList>() );
reduce( cellVolume, sumOp<scalarList>() );
reduce( volume, sumOp<scalarList>() );
reduce( position, sumOp<vectorList>() );
reduce( velocity, sumOp<vectorList>() );

//-ML: Print out Info
for(int i=1; i<=maxID; i++)
{
    Info << "/************** Structure Number is " << i << " **************/" << endl;
    Info << "Number of Cells: " << noCells[i] << endl;
    Info << "Volume of Structure: " << volume[i] << endl;
    Info << "position of Structure: " << position[i] / volume[i] << endl;
    Info << "velocity of Structure: " << velocity[i] / volume[i] << endl;
}


//-ML: Process small structures
//-ML: Add code to inject properties (e.g., position, velocity) into the cloud
for (label i = 1; i <= maxID; i++)
{
    if (noCells[i] < minCells)
    {
        
        forAll(mesh.cells(), cellI)
        {
            label strcutureID = floor(vofID_[cellI] + 0.5);
            if (strcutureID == i)
            {
                alpha1[cellI] = 1.0; //-ML: Set alpha1 to 1 for small structures
            }
        }
        
        //-ML: Calculate center position of the structure
        vector bubblePosition = position[i] / volume[i] ; 
        //-ML: Find the cell that contains this position
        label cellSample = mesh.findCell(bubblePosition); 
        //-ML: Calculate average velocity of the structure
        vector bubbleVelocity = velocity[i] / volume[i] ; 

        if (cellSample > -1)
        {
            //-ML: Apply corrections to position for 2-D cases
            meshTools::constrainToMeshCentre(mesh, bubblePosition);
            //-ML: Inject into Lagrangian cloud
            basicKinematicBubbleParcel* pPtr = new basicKinematicBubbleParcel(mesh, bubblePosition, cellSample); 
            //-ML: Number of particles per parcel
            pPtr->nParticle() = 1; 
            //-ML: Particle diameter
            pPtr->d() = Foam::cbrt(6.0 * volume[i] / constant::mathematical::pi); 
            //-ML: Velocity
            pPtr->U() = bubbleVelocity;
             
            basicKinematicBubbleParcel::trackingData td(parcels);
            scalar trackTime = mesh.time().deltaTValue();
            td.part() = basicKinematicBubbleParcel::trackingData::tpLinearTrack;

            //-ML:  Check/set parcel properties
            parcels.setParcelThermoProperties(*pPtr, trackTime);
            parcels.checkParcelProperties(*pPtr, trackTime, false);

            //-ML:  Apply correction to velocity for 2-D cases
            meshTools::constrainDirection
            (
                mesh,
                mesh.solutionD(),
                pPtr->U()
            );
            
            if (pPtr->move(parcels, td, trackTime)) //-ML: Check the possiblity of using move function
            {
                Pout << "Injecting Structure " << i 
                << " into Lagrangian cloud by processor: " << Pstream::myProcNo() 
                << " with velocity of " << pPtr->U() 
                << " and diamter " << pPtr->d() << endl;
                parcels.addParticle(pPtr); //-ML: Add particle to the cloud
            }
            else
            {
                delete pPtr;
            }
            
        }
        
    }
}

