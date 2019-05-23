% Script to convert from XML to .DATA (LAMMPS)
% For now it works with amorphous homopolymers
% 2019-05-23
%@author: Enzozr
clear all 
%% Open and read file

n1='HomopolymerRelaxed';
n2='.xml';

name_file=strcat(n1,n2);
XML_FILE=fopen(name_file);
FILE_OPEN=1;

while FILE_OPEN==1
    
    string=(strsplit(fgetl(XML_FILE)));
    
    %After line "<box.." there are the box dimensions
    if(strcmpi(string{1},'<box'));
        box.x=strsplit(string{2},'"');
        box.y=strsplit(string{2},'"');
        box.z=strsplit(string{2},'"');
        box.x=str2num(box.x{2});
        box.y=str2num(box.y{2});
        box.z=str2num(box.z{2});
    end
    
    %atoms positions
    if(strcmpi(string{1},'<position'));
        N_atoms=strsplit(string{end},'"');
        N_atoms=str2num(N_atoms{2});

        for i=1:1:N_atoms
            a=str2num(fgetl(XML_FILE));
            for k=1:1:3
                data_atoms(i,1)=i;
                data_atoms(i,3+k)=a(k);
            end
        end
    end
    
    %image number
    if(strcmpi(string{1},'<image'));
        for i=1:1:N_atoms
            a=str2num(fgetl(XML_FILE));
            for k=1:1:3
                data_atoms(i,6+k)=a(k);
            end
        end
    end
    
    %atoms velocities
    if(strcmpi(string{1},'<velocity'));
        for i=1:1:N_atoms
            a=str2num(fgetl(XML_FILE));
            for k=1:1:3
                data_velocities(i,1)=i;
                data_velocities(i,1+k)=a(k);
            end
        end
    end
    
    %type A->1
    %type B->2
    %another type ->99
    if(strcmpi(string{1},'<type'));
        for i=1:1:N_atoms
            string=(strsplit(fgetl(XML_FILE)));
            if(strcmpi(string,'A'))
                data_atoms(i,3)=1;
            elseif(strcmpi(string,'B'))
                data_atoms(i,3)=2;
            else 
                data_atoms(i,3)=99;
            end
        end
    end        
    
    %bonds information
    if(strcmpi(string{1},'<bond'));
        N_bonds=strsplit(string{end},'"');
        N_bonds=str2num(N_bonds{2});  
        for j=1:1:N_bonds
            string=(strsplit(fgetl(XML_FILE)));
            data_bonds(j,1)=str2num(string{2});
            data_bonds(j,2)=str2num(string{3});
        end
        %XML reading finished, close the file
        FILE_OPEN=0;
    end  
end
fclose all;

%In Hoomd atoms starts from 0
data_bonds=data_bonds+1;

clear a FILE_OPEN i j k name_file string XML_FILE
%% Add Number of molecule to atoms
N_molecule=1;
a1=data_bonds(1,1);
b1=data_bonds(1,2);
%The first atom...
data_atoms(1,2)=1;
%i is the atom number
i=1;
for j=2:1:N_bonds
    a2=data_bonds(j,1);
    b2=data_bonds(j,2);
    %when two consecutive atoms aren't bonded
    if ne(a2,b1)
        i=i+1;
        data_atoms(i,2)=N_molecule;
        %the next atom belongs to a new molecule
        N_molecule=N_molecule+1;
    end
    a1=a2;
    b1=b2;
    i=i+1;
    data_atoms(i,2)=N_molecule;
end
data_atoms(N_atoms,2)=N_molecule;
%% Clear
clear a1 a2 b1 b2 i j
save workspace.mat
%% Write .data
clear DATA_TEXT

%heading
DATA_TEXT{1,1}='LAMMPS data file via XML, version 31 Mar 2017, timestep = 0';
DATA_TEXT{end+1,1}=' ';
DATA_TEXT{end+1,1}=strcat(num2str(N_atoms),' atoms');
DATA_TEXT{end+1,1}=strcat(num2str(max(data_atoms(:,3))),' atom types');
DATA_TEXT{end+1,1}=strcat(num2str(N_bonds),' bonds');
DATA_TEXT{end+1,1}=strcat('1',' bond types');
DATA_TEXT{end+1,1}=' ';
%char(...) is to convert from cell to string
DATA_TEXT{end+1,1}=char(strcat(num2str(-box.x/2),{' '},num2str(box.x/2),' xlo xhi'));
DATA_TEXT{end+1,1}=char(strcat(num2str(-box.y/2),{' '},num2str(box.y/2),' ylo yhi'));
DATA_TEXT{end+1,1}=char(strcat(num2str(-box.z/2),{' '},num2str(box.z/2),' zlo zhi'));
DATA_TEXT{end+1,1}=' ';
DATA_TEXT{end+1,1}='Masses';
DATA_TEXT{end+1,1}=' ';
DATA_TEXT{end+1,1}='1 1';
DATA_TEXT{end+1,1}=' ';
DATA_TEXT{end+1,1}='Pair Coeffs # lj/cut/gpu';
DATA_TEXT{end+1,1}=' ';
DATA_TEXT{end+1,1}='1 1 1';
DATA_TEXT{end+1,1}=' ';
DATA_TEXT{end+1,1}='Bond Coeffs # fene';
DATA_TEXT{end+1,1}=' ';
DATA_TEXT{end+1,1}='1 30 1.5 1 1';
DATA_TEXT{end+1,1}=' ';

DATA_TEXT{end+1,1}='Atoms # molecular';
DATA_TEXT{end+1,1}=' ';
for i=1:1:N_atoms
    DATA_TEXT{end+1,1}=char(strcat(num2str(data_atoms(i,1)),{' '},num2str(data_atoms(i,2)),{' '},num2str(data_atoms(i,3)),{' '},num2str(data_atoms(i,4)),{' '},num2str(data_atoms(i,5)),{' '},num2str(data_atoms(i,6)),{' '},num2str(data_atoms(i,7)),{' '},num2str(data_atoms(i,8)),{' '},num2str(data_atoms(i,9))));
end
DATA_TEXT{end+1,1}=' ';
DATA_TEXT{end+1,1}='Velocities';
DATA_TEXT{end+1,1}=' ';
for i=1:1:N_atoms
    DATA_TEXT{end+1,1}=char(strcat(num2str(data_atoms(i,1)),{' '},num2str(data_velocities(i,2)),{' '},num2str(data_velocities(i,3)),{' '},num2str(data_velocities(i,4))));
end
DATA_TEXT{end+1,1}=' ';
DATA_TEXT{end+1,1}='Bonds';
DATA_TEXT{end+1,1}=' ';
for j=1:1:N_bonds
    DATA_TEXT{end+1,1}=char(strcat(num2str(j),{' '},'1',{' '},num2str(data_bonds(j,1)),{' '},num2str(data_bonds(j,2))));
end

lmp_file=fopen(strcat(n1,'.data'),'wt');
fprintf(lmp_file,'%s\n',DATA_TEXT{:});
fclose(lmp_file);

clear i j lmp_file