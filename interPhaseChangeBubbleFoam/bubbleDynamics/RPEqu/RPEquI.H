
inline Foam::RPEqu::RPEqu
(
)
:
    ODESystem(),
    pL_(0),
    rhoL_(0),
    muL_(0),
    sigmaL_(0),
    p0_(0),      
    R0_(0),
    pv_(0),
    kappa_(0)
{
}

inline Foam::RPEqu::~RPEqu()
{}


//-ML: Member functions
inline void Foam::RPEqu::setValues
(
    const scalar& pL,
    const scalar& rhoL,
    const scalar& muL,
    const scalar& sigmaL,
    const scalar& p0,      
    const scalar& R0,
    const scalar& pv,
    const scalar& kappa
)
{
    pL_ = pL;
    rhoL_=  rhoL;
    muL_= muL;
    sigmaL_ = sigmaL;
    p0_ = p0;      
    R0_ = R0;
    pv_ = pv;
    kappa_ = kappa;
}

inline void Foam::RPEqu::derivatives
(
    const scalar x,
    const scalarField& y,
    scalarField& dydx
) const
{
    scalar nuL=(muL_/rhoL_);
    scalar kappa1=kappa_;
    if (y[0] >= R0_) 
    {
        kappa1=1.0 ;
    }
    scalar pB=(p0_+(2*sigmaL_/R0_)-pv_)*(Foam::pow((R0_/y[0]),3*kappa1))+pv_;
    scalar dP=pL_-pB;

    dydx[0] = y[1];
    dydx[1] = -((3/2)*(sqr(y[1])/y[0]))-((4*nuL)*(y[1]/sqr(y[0])))-((2*sigmaL_)/(rhoL_*sqr(y[0])))-(dP/(rhoL_*y[0]));
}


inline void Foam::RPEqu::jacobian
(
    const scalar x,
    const scalarField& y,
    scalarField& dfdx,
    scalarSquareMatrix& dfdy
) const
{
    scalar nuL=(muL_/rhoL_);
    scalar kappa1=kappa_;
    if (y[0] >= R0_)
    {
        kappa1=1.0 ;
    }
    scalar pB=(p0_+(2*sigmaL_/R0_)-pv_)*(Foam::pow((R0_/y[0]),3*kappa1))+pv_;
    scalar dP=pL_-pB;

    dfdx[0] = 0.0;
    dfdx[1] = 0.0;
    dfdy[0][0] = 0.0;
    dfdy[0][1] = 1.0;
    dfdy[1][0] = ((3/2)*(sqr(y[1])/sqr(y[0])))+((8*nuL*y[1])/(pow(y[0],3)))+((4*sigmaL_)/(rhoL_*pow(y[0],3)))+(dP/(rhoL_*sqr(y[0])));
    dfdy[1][1] = -((3*y[1])/y[0])-((4*nuL)/sqr(y[0]));
}
