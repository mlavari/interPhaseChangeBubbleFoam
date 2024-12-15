//-ML: Set data from bubble
vector posB = this->position();
scalar radB = this->d() / 2;
label startCell = p.cell();
scalar totalBubbleVolume = (p.nParticle()*p.volume());

//-ML: DynamicList to store cells and their distances
DynamicList<std::pair<label, scalar>> cellsWithDistances;

//-ML: Add the start cell to the list 
vector cellCenter = this->mesh().cellCentres()[startCell];
cellsWithDistances.append(std::make_pair(startCell, mag(cellCenter - posB)));  //-ML: Add the first cell with distance

//-ML: Create a set to avoid re-visiting cells
labelHashSet visitedCells;
visitedCells.insert(startCell);

//-ML: Dynamic list to store all neighbors for further exploration
DynamicList<label> cells(0);
cells.append(this->mesh().cellCells()[startCell]);

//-ML: Start outer loop to explore neighbors
while (cells.size() > 0)
{
    //-ML: New dynamic list to store newly found neighbors within the radius
    DynamicList<label> newNeighbours;

    //-ML: Loop over all current neighbors in the cells list
    forAll(cells, i)
    {
        label cellI = cells[i];
        if (!visitedCells.found(cellI))  //-ML: Check if cellI is already visited
        {
            cellCenter = this->mesh().cellCentres()[cellI];

            //-ML: Check if this neighbor cell's center is within the radius
            scalar distance = mag(cellCenter - posB);  //-ML: Calculate distance to posB

            if (distance <= (2.0 * radB))
            {
                cellsWithDistances.append(std::make_pair(cellI, distance));  //-ML: Store the cell and its distance
                newNeighbours.append(this->mesh().cellCells()[cellI]);  //-ML: Add neighbors to the list
            }

            visitedCells.insert(cellI);  //-ML: Mark cell as visited
        }
    }

    //-ML: Sort the cells based on their distance to posB
    std::sort(cellsWithDistances.begin(), cellsWithDistances.end(),
              [](const std::pair<label, scalar>& a, const std::pair<label, scalar>& b)
              {
                  return a.second < b.second;  //-ML: Compare distances
              });

    //-ML: Continue to the next iteration with new Neighbours
    cells.clear();
    cells.transfer(newNeighbours);
}

//-ML: Initialize the boolean flag for proximity to the interface
bool nearInterface = false;

scalar alphaThreshold = this->LE_alphaThreshold(); //-ML: Set threshold for interface proximity


//-ML: Check if cellsWithDistances contains only one cell
if (cellsWithDistances.size() == 1)
{
    //-ML: Initialize a new DynamicList to hold the single cell and its neighbors
    DynamicList<label> singleCellAndNeighbors;

    //-ML: Get the single cell and add it to the new list
    label onlyCell = cellsWithDistances[0].first;
    singleCellAndNeighbors.append(onlyCell);

    //-ML: Add the neighbors of this single cell to the list
    singleCellAndNeighbors.append(this->mesh().cellCells()[onlyCell]);
 
    //-ML: Loop over the single cell and its neighbors to check alpha values
    forAll(singleCellAndNeighbors, j)
    {
        label cellJ = singleCellAndNeighbors[j];
        scalar alphaValue = cloud.alphal()[cellJ];  //-ML: Retrieve the alpha value for cellJ

        //-ML: Check if alpha value exceeds the threshold
        if (alphaValue < alphaThreshold)
        {
            nearInterface = true;
            break;  //-ML: Exit loop early if the condition is met
        }
    }
}
else
{
    //-ML: loop over cellsWithDistances and check alpha values
    forAll(cellsWithDistances, j)
    {
        label cellI = cellsWithDistances[j].first;  //-ML: Extract cell label
        scalar alphaValue = cloud.alphal()[cellI];  //-ML: Retrieve the alpha value for cellI

        //-ML: Check if alpha value exceeds threshold
        if (alphaValue < alphaThreshold)
        {
            nearInterface = true;
            break;  //-ML: Exit loop early if the condition is met
        }
    }
}

//-ML: Gaussian standard deviation
scalar sigma = radB / this->LE_deviation(); 

scalar totalWeightWithinRadius = 0.0;

scalar minCellOccupancy = this->LE_minCellOccupancy();

//-ML: Maximum fraction of volume for any single cell
scalar maxFraction = 0.0; 

//-ML: Accumulate excess volume to redistribute
scalar excessVolume = 0.0; 

//-ML: Set a threshold for the number of cells within the radius
label cellThreshold = this->LE_cellThreshold(); 

//-ML: Minimum bubble radius threshold 
scalar bubbleSizeThreshold = this->LE_bubbleSizeThreshold(); 

//-ML: Calculate the Gaussian weight for each cell and accumulate total weight within radius
List<scalar> gaussianWeights(cellsWithDistances.size());

forAll(cellsWithDistances, j) {

    scalar gaussFactor;

    if (cellsWithDistances.size() == 1) {
        gaussFactor = 1.0;  //-ML: Set gaussFactor to 1 if there is only one cell
    } else {
        scalar distance = cellsWithDistances[j].second;  //-ML: Get distance from pair
        //-ML: Calculate Gaussian weight for this cell
        gaussFactor = exp(-0.5 * pow(distance / sigma, 2.0));
    }
    gaussianWeights[j] = gaussFactor;

    //-ML: Accumulate total weight
    totalWeightWithinRadius += gaussFactor;
}

//-ML: Create a list to store cell IDs and distances with excess volume
DynamicList<label> excessVolumeCells;

//-ML: Distribute the bubble volume using normalized Gaussian weights
forAll(cellsWithDistances, j) {

    label cellJ = cellsWithDistances[j].first;
    //-ML:scalar distance = cellsWithDistances[j].second;  //-ML: Get distance from pair
    scalar cellVolume = this->mesh().cellVolumes()[cellJ];

    //-ML: Normalize the Gaussian weight and compute `vi_j` for this cell
    scalar normalizedWeight = gaussianWeights[j] / totalWeightWithinRadius;

    //-ML:Info << "normalizedWeight" << normalizedWeight << endl;
    scalar vi_j = normalizedWeight * totalBubbleVolume / cellVolume;

    maxFraction = (cloud.betal()[cellJ] - minCellOccupancy);
    //-ML: Apply the limit on maximum fraction of volume
    scalar cappedVi_j = min(vi_j, maxFraction);

    //-ML: Calculate any excess volume that needs redistribution
    if (vi_j > maxFraction) {
        excessVolume += (vi_j - maxFraction) * cellVolume;
    } else if (vi_j < maxFraction) {
        excessVolumeCells.append(cellJ);  //-ML: Store the cell 
    }

    //-ML: Update beta value with limit check
    cloud.betal()[cellJ] -= cappedVi_j;
    
    //-ML: Check if the updated beta value goes below 1; if so, set alpha to 1
    if (cloud.betal()[cellJ] < 1 && !nearInterface) 
    {
        cloud.alphal()[cellJ] = 1.0;
    }
}

//-ML: Check if cellsWithDistances has only one cell and
// there is still excessVolume, so there is no excessVolumeCells
if (cellsWithDistances.size() == 1 && excessVolume > 0) {

    label cellJ = cellsWithDistances[0].first;
    //-ML: Add the neighbors of this single cell to the list
    excessVolumeCells.append(this->mesh().cellCells()[cellJ]);    

}

//-ML: Redistribute any excess volume among cells that still have available capacity
forAll(excessVolumeCells, j) {

    label cellJ = excessVolumeCells[j];  //-ML: Get the cell ID

    scalar cellVolume = this->mesh().cellVolumes()[cellJ];

    //-ML: Log the initial beta value before redistribution
    scalar initialBeta = cloud.betal()[cellJ];

    //-ML: Determine remaining capacity for this cell
    scalar remainingCapacity = initialBeta - minCellOccupancy;

    //-ML: Only distribute to cells that still have capacity available
    if (remainingCapacity > 0 && excessVolume > 0) {
        //-ML: Calculate the portion of excess volume this cell can take
        scalar additionalVolume = min(excessVolume, remainingCapacity * cellVolume);

        //-ML: Apply additional volume to beta
        cloud.betal()[cellJ] -= additionalVolume / cellVolume;

        //-ML: Reduce excess volume by the amount distributed
        excessVolume -= additionalVolume;

        if (cloud.betal()[cellJ] < 1 && !nearInterface) 
        {
            cloud.alphal()[cellJ] = 1.0;
        }

    }
}

//-ML: Log the information
Info << "This bubble with ID: " << p.origId() 
<< "\n with the processor ID: " << p.origProc()
<< "\n in position: " << posB 
<< "\n with radius: " << radB
<< "\n within center cellID: " << startCell
<< "\n with bubble volume: " << totalBubbleVolume
<< "\n final excessVolume: " << excessVolume 
<< "\n the number of cells within 2xdiameter: " <<  cellsWithDistances.size() << endl;

//-ML: Apply the condition for number of cells within the radius or proximity to interface
if ((cellsWithDistances.size() > cellThreshold || nearInterface) && radB >= bubbleSizeThreshold) 
{
    forAll(cellsWithDistances, j) {
        label cellJ = cellsWithDistances[j].first;  //-ML: Get the cell ID
        //-ML: Decrease alphal by betal, ensuring alphal doesn't go below zero
        cloud.alphal()[cellJ] = max(cloud.alphal()[cellJ] - (1.0 - cloud.betal()[cellJ]), 0.0);
        //-ML: Set betal to 1
        cloud.betal()[cellJ] = 1.0;
    }

    //-ML: Stop tracing the bubble by setting it as inactive
    td.keepParticle = false;  
    Info << "This bubble tracing deactivated for particle at position: " << posB << endl;

    //-ML: Log the reason for deactivating bubble tracing
    if (cellsWithDistances.size() > cellThreshold) {
        Info << "\n Reason: Number of cells within the radius exceeded the threshold." << endl;
    } 
    else if (nearInterface) {
        Info << "\n Reason: Proximity to interface detected." << endl;
    }
}
else if (radB < bubbleSizeThreshold) 
{
    //-ML: If bubble size is below the threshold, deactivate tracking
    td.keepParticle = false;
    Info << "This bubble tracing deactivated for particle at position: " << posB << endl;
    Info << "\n Reason: Bubble radius (" << radB << ") is below the size threshold (" << bubbleSizeThreshold << ")." << endl;
}

//-ML: Activation control for box boundary check
bool boxCheckEnabled = this->LE_boxCheckEnabled(); 

//-ML: Define the box boundaries
vector boxTopLeftCorner = this->LE_boxTopLeftCorner();  //-ML: the top-left corner coordinates
vector boxBottomRightCorner = this->LE_boxBottomRightCorner();  //-ML: the bottom-right corner coordinates

//-ML: Check if the bubble parcel is within the box boundaries
bool withinBox = (posB.x() >= boxTopLeftCorner.x() && posB.x() <= boxBottomRightCorner.x() &&
                  posB.y() <= boxTopLeftCorner.y() && posB.y() >= boxBottomRightCorner.y() &&
                  posB.z() >= boxTopLeftCorner.z() && posB.z() <= boxBottomRightCorner.z());

//-ML: If bubble parcel is outside the box and box check is enabled, deactivate it
if (!withinBox && boxCheckEnabled) {
    td.keepParticle = false;  
    Info << "This bubble at position " << posB << " is outside the box and removed." << endl;
} 
