//-ML: Set tracking values from cells
this->setCellValues(cloud, td, this->averagingDistance());
//-ML: Create the selected ODE system solver
RPEqu_.setValues(
                td.pc(),
                td.rhoc(),
                td.muc(),
                this->bubbleSigma(),
                this->p0(),
                this->R0(),
                this->pv(),
                this->bubbleKappa()
                );


dictionary solverType;
solverType.add("solver", this->ODESolverType());

autoPtr<ODESolver> RPEquSolver = ODESolver::New(RPEqu_, solverType);
scalar xStart = this->age();
const scalar dx = this->RPdT();

//-ML: Initial radius and radiusGrowthRate
scalarField yStart(RPEqu_.nEqns());
yStart[0] = (this->d()/2);
yStart[1] = this->R_dot_;

//-ML: Integration initial step
scalar dxEst = 1e-10;
scalar xEnd = 0;

//-ML: Required to store dydx
scalarField dyStart(RPEqu_.nEqns());
const label n = dt/dx;

//-ML: Integration loop
for (label i=0; i<n; i++)
{
    xEnd = xStart + dx;
    RPEqu_.derivatives(xStart, yStart, dyStart);
    RPEquSolver->relTol() = 1e-5; //-ML:For runge kuta
    RPEquSolver->solve(xStart, xEnd, yStart, dxEst);
    xStart = xEnd;
}

this->d_ = (yStart[0]*2);
this->R_dot_ = yStart[1];

this->setCellValues(cloud, td);
