function [Mesh, U] = ImportOldEbbSolution(file)

	f = fopen(file, 'rb');
	
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
	
	nCells = fread(f, 1, 'uint64', 'b');
	U = fread(f, [4, nCells], 'double', 'b')';
	
	Mesh.V = nodes;
	Mesh.E2N = e2n;
	Mesh.IE = ie;
	Mesh.BE = be;
	Mesh.BName = tags;
	
end