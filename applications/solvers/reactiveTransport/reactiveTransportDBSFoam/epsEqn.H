{
    Info << "R=" << sum(R*mesh.V()) << endl;
    solve
    (
        fvm::ddt(eps)
        ==
        R
    );
    
   eps = min(max(1e-4,eps),1.0);
   Info << "fluid fraction" 
        << "  Min(eps) = " << gMin(eps)
        << "  Max(eps) = " << gMax(eps)
        << endl;


}
