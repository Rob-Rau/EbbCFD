function [t, forces] = ImportEbbForces(forceFile)

	f = fopen(forceFile, 'rb');
	
	forces = zeros(2,2);
	t = zeros(2,1);
	
	i = 1;
	while(0 == feof(f))
		tmp = fread(f, 1, 'double', 'b');
		if(0 == feof(f))
			t(i) = tmp;
			forces(i,:) = fread(f, 2, 'double', 'b');
		end
		i = i + 1;
	end
end