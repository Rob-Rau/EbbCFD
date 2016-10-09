function [Mesh, U, t, dt] = ImportEbbSolution(meshFile, slnFile)

	f = fopen(meshFile, 'rb');
	
	nNodes = fread(f, 1, 'uint64', 'b');
	nodes = fread(f, [2, nNodes], 'double', 'b')';
	
	nEls = fread(f, 1, 'uint64', 'b');
	e2n = fread(f, [3, nEls], 'double', 'b')';
	
	ieSize = fread(f, 1, 'uint64', 'b');
	ie = fread(f, [4, ieSize], 'double', 'b')';
	
	beSize = fread(f, 1, 'uint64', 'b');
	be = fread(f, [4, beSize], 'double', 'b')';
	
	nTags = fread(f, 1, 'uint64', 'b');
	for i = 1:nTags
		tagLen = fread(f, 1, 'uint32', 'b');
		tags{i} = fread(f, tagLen, '*char');
	end
	fclose(f);
	
	f = fopen(slnFile, 'rb');
	
	slnMagic = fread(f, 1, 'uint32', 'b');
	slnVer = fread(f, 1, 'uint32', 'b');
	nCells = fread(f, 1, 'uint32', 'b');
	t = fread(f, 1, 'double', 'b');
	dt = fread(f, 1, 'double', 'b');
	
	U = fread(f, [4, nCells], 'double', 'b')';
	
	fclose(f);
	
	Mesh.V = nodes;
	Mesh.E2N = e2n;
	Mesh.IE = ie;
	Mesh.BE = be;
	Mesh.BName = tags;
	
end