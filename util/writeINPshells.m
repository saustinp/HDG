function writeINPshells(fileName, pos, elem, elem2sol, elemfield)
%WRITEINPSHELLS : Write a INP file format for Paraview
%   Detailed explanation goes here


% Number of points, dimension elements and vertex per element
[np,nd]  = size(pos);
[ne,npv] = size(elem);

% Get element type
if npv==3
    elementTypeString = 'tri';
elseif npv==4 && nd==2
    elementTypeString = 'quad';
elseif npv==4 && nd==3
    elementTypeString = 'tet';
elseif npv==6 && nd==3
    elementTypeString = 'prism';
elseif npv==8
	elementTypeString = 'hex';
else
    error('Type of element not available.');
end


file = fopen([fileName,'.inp'],'w');

% Write header
toWrite = num2str([np, ne, 0, ne, 0]); 
writeStringMatrix(file, toWrite);

% Save nodes Coordinates
toWrite = num2str([(1:np)', pos]);
writeStringMatrix(file, toWrite);

% Save Element Connectivities
toWrite = [num2str((1:ne)'), repmat('   ',[ne,1]), num2str(elem2sol), repmat(['   ' elementTypeString '   '],[ne,1]), num2str(elem)];
writeStringMatrix(file, toWrite);

% Writing the quads thickness
toWrite = ['1   ', num2str(1)];
writeStringMatrix(file, toWrite);
toWrite = 'Thickness, None';        % None refers to the units of the variable
writeStringMatrix(file, toWrite);
toWrite = num2str([(1:ne)', elemfield]);
writeStringMatrix(file, toWrite);

fclose(file);

end


function writeStringMatrix(file, toWrite)

numLines = size(toWrite,1);

for i=1:numLines
    str = toWrite(i,:);
    fprintf(file, '%s\n', str);
end

end


