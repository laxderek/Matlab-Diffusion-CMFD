classdef Mesh

    properties
    phi;
    dtilde;
    Leakage;
    Jsurf;
	X;
	Y;
	Z;
	I;
	J;
	K;
	INDEX;
	MAT;
	dxyz;
    D;
    sigS;
    sigA;
    sigT;
    nusigF;
    chi;
	albedo;
	mats;
	x;
    y;
    z;
    g;
    gridreduce;
    dhat;
    end
	
    methods
        
        %constructor
        function self=Mesh(x,y,z,g) 
            self.phi = zeros(x*y*z*g,1);
            self.Leakage = zeros(x,y,z,g,6);
            self.Jsurf = zeros(x,y,z,g,6);
            self.dtilde = zeros(x,y,z,g,6);
            self.dhat = zeros(x,y,z,g,6);
            self.X = zeros(x,y,z);
            self.Y = zeros(x,y,z);
            self.Z = zeros(x,y,z);
            self.I = zeros(x,y,z);
            self.J = zeros(x,y,z);
            self.K = zeros(x,y,z);
            self.D = zeros(x,y,z,g);
            self.sigA = zeros(x,y,z,g);
            self.sigS = zeros(x,y,z,g,g);
            self.sigT = zeros(x,y,z,g);
            self.nusigF = zeros(x,y,z,g);
            self.chi = zeros(x,y,z,g);
            self.INDEX = zeros(x,y,z);
            self.MAT = zeros(x,y,z);
            self.dxyz = zeros(x,y,z,3);
            self.gridreduce = zeros(3,1);
            for i = 1:x
                for j = 1:y
                    for k = 1:z
                        self.I(i,j,k) = i;
                        self.J(i,j,k) = j;
                        self.K(i,j,k) = k;
                    end
                end
            end			
            self.x = x;
            self.y = y;
            self.z = z;
            self.g = g;
            self.albedo = zeros(3,1);
        end
	
        function self=setNumMats(self,num)
            self.mats = zeros(num,1);%new Material[num];
        end
	
        function self=setMatAtLoc(self,i,j,k,mat,material) 

            self.MAT(i,j,k) = mat;        
            self.D(i,j,k,:) = material.D;
            self.sigT(i,j,k,:) = material.sigT;
            self.sigA(i,j,k,:) = material.sigA;
            self.sigS(i,j,k,:,:) = material.sigS;
            self.chi(i,j,k,:) = material.chi;
            self.nusigF(i,j,k,:) = material.nusigF;
            
        end
        
        function self=setAllMat(self,mat,material)
            for i = 1:self.x
                for j = 1:self.y
                    for k = 1:self.z
                        self.MAT(i,j,k) = mat;        
                        self.D(i,j,k,:) = material.D;
                        self.sigT(i,j,k,:) = material.sigT;
                        self.sigA(i,j,k,:) = material.sigA;
                        self.sigS(i,j,k,:,:) = material.sigS;
                        self.nusigF(i,j,k,:) = material.nusigF;
                        self.chi(i,j,k,:) = material.chi;
                    end
                end
            end
        end	
	
        function self=setBC(self,xBC, yBC, zBC)
            self.albedo(1) = xBC;
            self.albedo(2) = yBC;
            self.albedo(3) = zBC;
        end
	
        function self=setEqualDist(self,dim, lenx, leny, lenz)
            dx = lenx / self.x;
            totx = 0;
            for i = 1:self.x
                for j = 1:self.y
                    for k = 1:self.z
                        self.dxyz(i,j,k,1) = dx;
                        self.X(i,j,k) = totx;
                    end
                end
                totx = totx + dx;
            end
            dy = leny / self.y;
            toty = 0;
            for j = 1:self.y
                for i = 1:self.x
                    for k = 1:self.z
                        if (dim > 1) %//y
                            self.dxyz(i,j,k,2) = dy;
                            self.Y(i,j,k) = toty;
                        else
                            self.dxyz(i,j,k,2) = 1;
                            self.Y(i,j,k) = 1;
                        end
                    end
                end
                toty = toty + dy;
            end
            dz = lenz / self.z;
            totz = 0;
            for k = 1:self.z
                for j = 1:self.y
                    for i = 1:self.x
                        if (dim > 2) %//z
                            self.dxyz(i,j,k,3) = dz;
                            self.Y(i,j,k) = totz;
                        else
                            self.dxyz(i,j,k,3) = 1;
                            self.Z(i,j,k) = 1;
                        end
                    end
                end
                totz = totz + dz;
            end
        end
    end
    
    methods(Static)
        
        function [coarse] = fineToCoarse(fine, phi, scale, dim)

            coarse = 0;
            ng = fine.g;
            %check that matrix can be scaled
            if (mod(fine.x,scale(1)) == 0 && mod(fine.y,scale(2)) == 0 && mod(fine.z,scale(3)) == 0)
                
                coarse = Mesh(fine.x / scale(1), fine.y / scale(2), fine.z / scale(3), fine.g);
                coarse.gridreduce = scale;
                coarse.albedo = fine.albedo;
                info = [dim ng fine.x fine.y fine.z];
                info2 = [dim ng coarse.x coarse.y coarse.z];
                %loop over coarse mesh condensing each dimension
                %(scalex,scaley,scalez)
                for i = 1:coarse.x
                    for j = 1:coarse.y
                        for k = 1:coarse.z
                            
                            %set all counters to zero
                            totphi = zeros(ng,1);
                            DiffP = zeros(ng,1);
                            sigmaTP = zeros(ng,1);
                            sigmaAP = zeros(ng,1);
                            sigmaSP = zeros(ng);
                            nusigmaFP = zeros(ng,1);
                            chiP = zeros(ng,1);
                            volWeightedPhi = zeros(ng);
                            coarse.X(i,j,k) = fine.X(i,j,k)*scale(1);  
                            coarse.Y(i,j,k) = fine.Y(i,j,k)*scale(2);  
                            coarse.Z(i,j,k) = fine.Z(i,j,k)*scale(3);
                            coarse.dxyz(i,j,k,1) = fine.dxyz(i,j,k,1) * scale(1);
                            coarse.dxyz(i,j,k,2) = fine.dxyz(i,j,k,2) * scale(2);
                            coarse.dxyz(i,j,k,3) = fine.dxyz(i,j,k,3) * scale(3);
                            
                            %set coarse mesh cell currents 
                            %%%ONLY WORKS IN 1D
                            for g = 1:fine.g
                                coarse.Jsurf(i,j,k,g,1) = fine.Jsurf((i-1) * scale(1) + 1,1,1,1);
                                coarse.Jsurf(i,j,k,g,2) = fine.Jsurf((i) * scale(1),1,1,2);
                            end
                            
                            for finex = (i-1) * scale(1) + 1: (i) * scale(1)
                                for finey = j * scale(2): (k+1) * scale(2) - 1
                                    for finez = k * scale(3): (k+1) * scale(3) - 1
                                        if (finex > fine.x || finey > fine.y || finez > fine.z) 
                                            continue;
                                        end
                                        
                                        
                                        for group  = 1:fine.g
                                        
                                            volWeightedPhi(group) = volWeightedPhi(group) + phi(indexToMat(finex,finey,finez,group,info)) * prod(fine.dxyz(i,j,k,:));                 
                                            totphi(group) = totphi(group) + phi(indexToMat(finex,finey,finez,group,info));
                                            DiffP(group) = DiffP(group) + fine.D(finex,finey,finez,group) * phi(indexToMat(finex,finey,finez,group,info));
                                            sigmaTP(group) = sigmaTP(group) + fine.sigT(finex,finey,finez,group) * phi(indexToMat(finex,finey,finez,group,info));
                                            sigmaAP(group) = sigmaAP(group) + fine.sigA(finex,finey,finez,group) * phi(indexToMat(finex,finey,finez,group,info));
                                            nusigmaFP(group) = nusigmaFP(group) + fine.nusigF(finex,finey,finez,group) * phi(indexToMat(finex,finey,finez,group,info));                                                
                                            chiP(group) = chiP(group) + fine.chi(finex,finey,finez,group) * phi(indexToMat(finex,finey,finez,group,info));                                                

                                            for h = 1:fine.g
                                                sigmaSP(group,h) = sigmaSP(group,h) + fine.sigS(finex,finey,finez,group,h) * phi(indexToMat(finex,finey,finez,group,info));
                                            end
                                            
                                        end
                                    end
                                end
                            end
                            for group  = 1:fine.g
                                    coarse.phi(indexToMat(i,j,k,group,info2)) = volWeightedPhi(group) / prod(coarse.dxyz(i,j,k,:));
                                    coarse.D(i,j,k,group) = DiffP(group) ./  totphi(group);
                                    coarse.sigT(i,j,k,group) = sigmaTP(group) ./ totphi(group);
                                    coarse.sigA(i,j,k,group) = sigmaAP(group) ./ totphi(group);
                                    coarse.nusigF(i,j,k,group) = nusigmaFP(group) ./ totphi(group);
                                    coarse.chi(i,j,k,group) = chiP(group) ./ totphi(group);
                                for h = 1:fine.g
                                    coarse.sigS(i,j,k,group,h) = sigmaSP(group,h) ./ totphi(group);
                                end
                            end
                        end
                    end
                end
            else
                print "Could not scale mesh properly"                
            end
        end
		
		
		
		
    end

end
