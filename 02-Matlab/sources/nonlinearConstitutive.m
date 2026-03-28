function [Cmat,Smat,W0] = nonlinearConstitutive(param,Fmat)
	% neo-Hookean material
	% 

	% right Cauchy-Green tensor
	C = Fmat'*Fmat;
	invC = inv(C);

	I1 = trace(C); I3 = det(C);

	lmbda = param(2) -  2/3*param(1);
	mu = param(1) - lmbda/2*log(I3);

	Iden = speye(3);

	% PK2 stress
	S = lmbda/2*log(I3)*invC + param(1)*(Iden - invC);
	Smat = zeros(4,4);
	Smat(1:2,1:2) = S(1:2,1:2);
	Smat(3:end,3:end) = S(1:2,1:2);

	% 4th order constitutive tensor
	Cmat = zeros(3,3);
	Cmat(1,1) = lmbda*invC(1,1)*invC(1,1) + mu*(invC(1,1)*invC(1,1) + invC(1,1)*invC(1,1));
	Cmat(1,2) = lmbda*invC(1,1)*invC(2,2) + mu*(invC(1,2)*invC(1,2) + invC(1,2)*invC(1,2));
	Cmat(1,3) = lmbda*invC(1,1)*invC(1,2) + mu*(invC(1,1)*invC(1,2) + invC(1,2)*invC(1,1));
	% 
	Cmat(2,1) = Cmat(1,2);
	Cmat(2,2) = lmbda*invC(2,2)*invC(2,2) + mu*(invC(2,2)*invC(2,2) + invC(2,2)*invC(2,2));
	Cmat(2,3) = lmbda*invC(2,2)*invC(1,2) + mu*(invC(2,1)*invC(2,2) + invC(2,2)*invC(2,1));
	% 
	Cmat(3,1) = Cmat(1,3);
	Cmat(3,2) = Cmat(2,3);
	Cmat(3,3) = lmbda*invC(1,2)*invC(1,2) + mu*(invC(1,1)*invC(2,2) + invC(1,2)*invC(2,1));

	% strain energy
	W0 = lmbda/8*(log(I3))^2 - param(1)/2*log(I3) + param(1)/2*(I1 - 3);

end
