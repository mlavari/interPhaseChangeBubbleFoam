#ifndef RPEqu_H
#define RPEqu_H

#include "argList.H"
#include "IOmanip.H"
#include "ODESystem.H"
#include "ODESolver.H"

namespace Foam
{

class RPEqu
:
    public ODESystem
{
private:
    // Private data members
    scalar pL_;
    scalar rhoL_;
    scalar muL_;    
    scalar sigmaL_;
    scalar p0_;      
    scalar R0_;
    scalar pv_;
    scalar kappa_;


public:

    // Constructor
    inline RPEqu();
    // Destructor
    inline ~RPEqu();
    

    // Member functions
    //  Set environment data  
    inline void setValues
    (
        const scalar& pL,
        const scalar& rhoL,
        const scalar& muL,
        const scalar& sigmaL,
        const scalar& p0,    
        const scalar& R0,
        const scalar& pv,
        const scalar& kappa
    );

    // Return a reference for dictionary
    inline label nEqns() const
    {
        return 2;
    }

    inline void derivatives
    (
        const scalar x,
        const scalarField& y,
        scalarField& dydx
    ) const;

    inline void jacobian
    (
        const scalar x,
        const scalarField& y,
        scalarField& dfdx,
        scalarSquareMatrix& dfdy
    ) const;    
};

};

#include "RPEquI.H"

#endif
