classdef Cell
	
    properties
	x;
    y;
    z;
	i;
    j;
    k;
    index; %x,y,z location in grid, i = 1D grid ID
	matid; %id of material in this cell
	dxyz; %lengths of dimensions of this cell
    end
    
    methods
     %   function self=Cell() 
     %       self.dxyz = zeros(3,1);
     %   end

        %function self=Cell(i, j, k, dxyz, index, matid)
    	function self=Cell()
            %self.i = i;
            %self.j = j;
            %self.k = k;
            %self.dxyz = dxyz;
            %self.index = index;
            %self.matid = matid;
            self.dxyz = zeros(3,1);
    	
        end

        function V = getVol(i)
            V = self.dxyz(1) * self.dxyz(2) * self.dxyz(3);
        end
    end
end
