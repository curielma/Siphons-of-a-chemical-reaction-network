--GOAL: computing siphons from a weakly reversible network
--INPUT: conservation law matrix A
--NOTES: network is specified after a choice of basis of ker(A);
--    	  constructed network is weakly connected
--CREDIT: Siphons in Chemical Reaction Networks
--     	   by Anne Shiu & Bernd Sturmfels
loadPackage "Polyhedra"

A = matrix{{0,0,1,1,1},{1,2,0,1,2}}

Q = coneFromVData A
F = facets Q

maxFaces = unique for i from 0 to numColumns F - 1 list(
    cols := {};
    j = 0; while j < numColumns A do(
    	if not contains(coneFromVData matrix F_i,coneFromVData matrix A_j) then cols = append(cols,j); 
	j = j+1;
	);
    unique flatten cols 
    )

R = QQ[x_1..x_(numColumns A)];
B = intersect(apply(maxFaces, i -> ideal apply(i, j -> x_(j+1))))

N = gens ker A
J = ideal for i from 0 to numColumns N - 1 list(
    gammaPlus = apply(entries N_i, j -> max(0,j));
    gammaMinus = - apply(entries N_i, j -> min(0,j));
    (product apply(#gammaPlus, j -> x_(j+1)^(gammaPlus_j))) - (product apply(#gammaMinus, j -> x_(j+1)^(gammaMinus_j)))
    )

saturate(J,B)
decompose oo
