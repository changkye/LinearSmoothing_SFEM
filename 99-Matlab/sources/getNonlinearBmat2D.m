function [Bmat,Bgeo,Fmat] = getNonlinearBmat2D(dNdX,wkU,nndof)
	% 

	% deformation gradient
	dxdX = wkU'*dNdX;
	Fmat = dxdX + eye(2);
	Fmat(3,3) = 1.0;

	% smoothed strain-displacement matrix
	Bmat = zeros(3,2*nndof);
	Bmat(1,1:2:end) = dNdX(:,1)*Fmat(1,1);
	Bmat(1,2:2:end) = dNdX(:,1)*Fmat(2,1);
	Bmat(2,1:2:end) = dNdX(:,2)*Fmat(1,2);
	Bmat(2,2:2:end) = dNdX(:,2)*Fmat(2,2);
	Bmat(3,1:2:end) = dNdX(:,2)*Fmat(1,1) + dNdX(:,1)*Fmat(1,2);
	Bmat(3,2:2:end) = dNdX(:,1)*Fmat(2,2) + dNdX(:,2)*Fmat(2,1);
	% 
	Bgeo = zeros(4,2*nndof);
	Bgeo(1,1:2:end) = dNdX(:,1);
	Bgeo(2,1:2:end) = dNdX(:,2);
	Bgeo(3,2:2:end) = dNdX(:,1);
	Bgeo(4,2:2:end) = dNdX(:,2);

end
