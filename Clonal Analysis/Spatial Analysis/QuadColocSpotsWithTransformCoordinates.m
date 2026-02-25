% Original script written by Mariano Gabitto. This version has been modified by Sophie Gobeil.

% <CustomTools>
%   <Menu>
%     <Item name="QuadColocSpotsWithTransformCoordinates" icon="Matlab" tooltip="triple coloc analysis">
%       <Command>MatlabXT::QuadColocSpotsWithTransformCoordinates(%i)</Command>
%     </Item>
%     
%   </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSpots">
%          <Item name="TripColocTransformCoords" icon="Matlab" tooltip="triple coloc analysis">
%            <Command>MatlabXT::QuadColocSpotsWithTransformCoordinates(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
% 
% </CustomTools>







function QuadColocSpotsWithTransformCoordinates( aImarisApplicationID )


global vImarisApplication;

if ~isa(aImarisApplicationID, 'Imaris.IApplicationPrxHelper')
  javaaddpath ImarisLib.jar
  vImarisLib = ImarisLib;
  if ischar(aImarisApplicationID)
    aImarisApplicationID = round(str2double(aImarisApplicationID));
  end
   vImarisApplication = vImarisLib.GetApplication(aImarisApplicationID);
else
  vImarisApplication = aImarisApplicationID;
end

vDataSet = vImarisApplication.GetDataSet.Clone;
vSurpassScene = vImarisApplication.GetSurpassScene;
aType = Imaris.tType.eTypeFloat;
% 
if strcmp(vDataSet.GetType,'eTypeUInt8')
    vImarisApplication.GetDataSet.SetType(aType);
elseif strcmp(vDataSet.GetType,'eTypeUInt16')
    vImarisApplication.GetDataSet.SetType(aType);
end    

% the user has to create a scene with some spots

if isequal(vSurpassScene, [])
    msgbox('Please create some Spots in the Surpass scene!');
    return;
end

prompt= {'Enter the maximum distance in Microns(um) for coloc','new Spot Radius','Counterclockwise rotation angle (degrees)'};
dlg_title='what a lovely day for counting spots';
num_lines=1;
def={'3','3', '0'};
UserInput=inputdlg(prompt,dlg_title,num_lines,def);
nSpotsObj=0;
MaxDist=str2double(UserInput{1})^2;
SpotRadius=str2double(UserInput{2});
Rotation=str2double(UserInput{3});
SpotRadius=round(SpotRadius);
for vChildIndex = 1:vSurpassScene.GetNumberOfChildren
    vDataItem = vSurpassScene.GetChild(vChildIndex - 1);
    if vImarisApplication.GetFactory.IsSpots(vDataItem)
        nSpotsObj=nSpotsObj+1;
        vSpots=vImarisApplication.GetFactory.ToSpots(vDataItem);
        vSpotValues{nSpotsObj}=struct2cell(struct(vSpots.Get));
        [vSpotsXYZ{nSpotsObj},vSpotsTime{nSpotsObj},vSpotsRadius{nSpotsObj}]=vSpotValues{nSpotsObj}{:};
        spotName{nSpotsObj}=char(vDataItem.GetName);
        spotColor{nSpotsObj}=vDataItem.GetColorRGBA;
        [Nspots(nSpotsObj),dim]=size(vSpotsXYZ{nSpotsObj});
    end
end  
display(spotName)
if nSpotsObj<4
    error('you must create 4 spots objects');
end   
if nSpotsObj>4
    if Nspots(5)==Nspots(6) && Nspots(6)==Nspots(7)
        if Nspots(5)==1
            goodStuff=true;
            ccPoint=vSpotsXYZ{5};
            bPoint =vSpotsXYZ{6};
            lPoint =vSpotsXYZ{7};
            tPoint =vSpotsXYZ{8};
        else
            error('number of spots in spot objects 5, 6, and 7 must be only 1');
        end
    else
        error('number of spots in spot objects 5, 6, and 7 must be only 1');
    end
else
    thisMessage='there are not enough spot objects to give transform coordinates.  These coordinate will be excluded from csv file';
    choice1=menu(thisMessage,'Continue','Cancel');
    if choice1==1
        excludeTransformPoints=true;
    else if choice1==2
            error('canceling extension')
        end
    end    
end    

% choice=menu('choose coloc result type', 'Position based coloc', 'object based coloc');
choice=2;
%Find coloc for spot1 and spot2
n1Colocs2=0;
loc1Co2=[];
Spot1Co2Logic=zeros(Nspots(1),1);
Spot2Co1Logic=zeros(Nspots(2),1);
if ~isempty(vSpotsXYZ{1}) || ~isempty(vSpotsXYZ{2})
    for i=1:Nspots(1)
        for j=1:Nspots(2)
            Dx=vSpotsXYZ{1}(i,1)-vSpotsXYZ{2}(j,1);
            Dy=vSpotsXYZ{1}(i,2)-vSpotsXYZ{2}(j,2);
            Dz=vSpotsXYZ{1}(i,3)-vSpotsXYZ{2}(j,3);
            Dist12=(Dx^2)+(Dy^2)+(Dz^2);
            if Dist12<=MaxDist
                Spot1Co2Logic(i)=1;
                Spot2Co1Logic(j)=1;
                n1Colocs2=n1Colocs2+1;
                loc1Co2(n1Colocs2,1)=mean([vSpotsXYZ{1}(i,1),vSpotsXYZ{2}(j,1)]);
                loc1Co2(n1Colocs2,2)=mean([vSpotsXYZ{1}(i,2),vSpotsXYZ{2}(j,2)]);
                loc1Co2(n1Colocs2,3)=mean([vSpotsXYZ{1}(i,3),vSpotsXYZ{2}(j,3)]);
                Stime12(n1Colocs2)=0;
                radii12(n1Colocs2)=SpotRadius;
            else
                if ~isequal(Spot1Co2Logic(i),1)
                    Spot1Co2Logic(i)=0;
                end
                
                if ~isequal(Spot2Co1Logic(j),1)
                    Spot2Co1Logic(j)=0;
                end
                 
            end    
        end
    end
end       
    
%Find coloc for spot1 and spot3
n1Colocs3=0;
loc1Co3=[];
Spot1Co3Logic=zeros(Nspots(1),1);
Spot3Co1Logic=zeros(Nspots(3),1);
if ~isempty(vSpotsXYZ{1}) || ~isempty(vSpotsXYZ{3})
    for i=1:Nspots(1)
        for j=1:Nspots(3)
            Dx=vSpotsXYZ{1}(i,1)-vSpotsXYZ{3}(j,1);
            Dy=vSpotsXYZ{1}(i,2)-vSpotsXYZ{3}(j,2);
            Dz=vSpotsXYZ{1}(i,3)-vSpotsXYZ{3}(j,3);
            Dist13=(Dx^2)+(Dy^2)+(Dz^2);
            if Dist13<=MaxDist
                Spot1Co3Logic(i)=1;
                Spot3Co1Logic(j)=1;
                n1Colocs3=n1Colocs3+1;
                loc1Co3(n1Colocs3,1)=mean([vSpotsXYZ{1}(i,1),vSpotsXYZ{3}(j,1)]);
                loc1Co3(n1Colocs3,2)=mean([vSpotsXYZ{1}(i,2),vSpotsXYZ{3}(j,2)]);
                loc1Co3(n1Colocs3,3)=mean([vSpotsXYZ{1}(i,3),vSpotsXYZ{3}(j,3)]);
                Stime13(n1Colocs3)=0;
                radii13(n1Colocs3)=SpotRadius;
            else
                if ~isequal(Spot1Co3Logic(i),1)
                    Spot1Co3Logic(i)=0;
                end    
            
                if ~isequal(Spot3Co1Logic(j),1)
                    Spot3Co1Logic(j)=0;
                end
     
            end    
        end
    end
end    

%Find coloc for spot1 and spot4
n1Colocs4=0;
loc1Co4=[];
Spot1Co4Logic=zeros(Nspots(1),1);
Spot4Co1Logic=zeros(Nspots(4),1);
if ~isempty(vSpotsXYZ{1}) || ~isempty(vSpotsXYZ{4})
    for i=1:Nspots(1)
        for j=1:Nspots(4)
            Dx=vSpotsXYZ{1}(i,1)-vSpotsXYZ{4}(j,1);
            Dy=vSpotsXYZ{1}(i,2)-vSpotsXYZ{4}(j,2);
            Dz=vSpotsXYZ{1}(i,3)-vSpotsXYZ{4}(j,3);
            Dist14=(Dx^2)+(Dy^2)+(Dz^2);
            if Dist14<=MaxDist
                Spot1Co4Logic(i)=1;
                Spot4Co1Logic(j)=1;
                n1Colocs4=n1Colocs4+1;
                loc1Co4(n1Colocs4,1)=mean([vSpotsXYZ{1}(i,1),vSpotsXYZ{4}(j,1)]);
                loc1Co4(n1Colocs4,2)=mean([vSpotsXYZ{1}(i,2),vSpotsXYZ{4}(j,2)]);
                loc1Co4(n1Colocs4,3)=mean([vSpotsXYZ{1}(i,3),vSpotsXYZ{4}(j,3)]);
                Stime14(n1Colocs4)=0;
                radii14(n1Colocs4)=SpotRadius;
            else
                if ~isequal(Spot1Co4Logic(i),1)
                    Spot1Co4Logic(i)=0;
                end
                
                if ~isequal(Spot4Co1Logic(j),1)
                    Spot4Co1Logic(j)=0;
                end
                 
            end    
        end
    end
end     

%Find coloc for spot2 and spot3
n2Colocs3=0;
loc2Co3=[];
Spot2Co3Logic=zeros(Nspots(2),1);
Spot3Co2Logic=zeros(Nspots(3),1);
if ~isempty(vSpotsXYZ{2}) || ~isempty(vSpotsXYZ{3})
    for i=1:Nspots(2)
        for j=1:Nspots(3)
            Dx=vSpotsXYZ{2}(i,1)-vSpotsXYZ{3}(j,1);
            Dy=vSpotsXYZ{2}(i,2)-vSpotsXYZ{3}(j,2);
            Dz=vSpotsXYZ{2}(i,3)-vSpotsXYZ{3}(j,3);
            Dist23=(Dx^2)+(Dy^2)+(Dz^2);
            if Dist23<=MaxDist
                Spot2Co3Logic(i)=1;
                Spot3Co2Logic(j)=1;
                n2Colocs3=n2Colocs3+1;
                loc2Co3(n2Colocs3,1)=mean([vSpotsXYZ{2}(i,1),vSpotsXYZ{3}(j,1)]);
                loc2Co3(n2Colocs3,2)=mean([vSpotsXYZ{2}(i,2),vSpotsXYZ{3}(j,2)]);
                loc2Co3(n2Colocs3,3)=mean([vSpotsXYZ{2}(i,3),vSpotsXYZ{3}(j,3)]);
                Stime23(n2Colocs3)=0;
                radii23(n2Colocs3)=SpotRadius;   
            else
                if ~isequal(Spot2Co3Logic(i),1)
                    Spot2Co3Logic(i)=0;
                end
              
                if ~isequal(Spot3Co2Logic(j),1)
                    Spot3Co2Logic(j)=0;
                end
               
            end    
        end
    end
end

%Find coloc for spot2 and spot4
n2Colocs4=0;
loc2Co4=[];
Spot2Co4Logic=zeros(Nspots(2),1);
Spot4Co2Logic=zeros(Nspots(4),1);
if ~isempty(vSpotsXYZ{2}) || ~isempty(vSpotsXYZ{4})
    for i=1:Nspots(2)
        for j=1:Nspots(4)
            Dx=vSpotsXYZ{2}(i,1)-vSpotsXYZ{4}(j,1);
            Dy=vSpotsXYZ{2}(i,2)-vSpotsXYZ{4}(j,2);
            Dz=vSpotsXYZ{2}(i,3)-vSpotsXYZ{4}(j,3);
            Dist24=(Dx^2)+(Dy^2)+(Dz^2);
            if Dist24<=MaxDist
                Spot2Co4Logic(i)=1;
                Spot4Co2Logic(j)=1;
                n2Colocs4=n2Colocs4+1;
                loc2Co4(n2Colocs4,1)=mean([vSpotsXYZ{2}(i,1),vSpotsXYZ{4}(j,1)]);
                loc2Co4(n2Colocs4,2)=mean([vSpotsXYZ{2}(i,2),vSpotsXYZ{4}(j,2)]);
                loc2Co4(n2Colocs4,3)=mean([vSpotsXYZ{2}(i,3),vSpotsXYZ{4}(j,3)]);
                Stime24(n2Colocs4)=0;
                radii24(n2Colocs4)=SpotRadius;   
            else
                if ~isequal(Spot2Co4Logic(i),1)
                    Spot2Co4Logic(i)=0;
                end
              
                if ~isequal(Spot4Co2Logic(j),1)
                    Spot4Co2Logic(j)=0;
                end
               
            end    
        end
    end
end

%Find coloc for spot3 and spot4
n3Colocs4=0;
loc3Co4=[];
Spot3Co4Logic=zeros(Nspots(3),1);
Spot4Co3Logic=zeros(Nspots(4),1);
if ~isempty(vSpotsXYZ{3}) || ~isempty(vSpotsXYZ{4})
    for i=1:Nspots(3)
        for j=1:Nspots(4)
            Dx=vSpotsXYZ{3}(i,1)-vSpotsXYZ{4}(j,1);
            Dy=vSpotsXYZ{3}(i,2)-vSpotsXYZ{4}(j,2);
            Dz=vSpotsXYZ{3}(i,3)-vSpotsXYZ{4}(j,3);
            Dist34=(Dx^2)+(Dy^2)+(Dz^2);
            if Dist34<=MaxDist
                Spot3Co4Logic(i)=1;
                Spot4Co3Logic(j)=1;
                n3Colocs4=n3Colocs4+1;
                loc3Co4(n3Colocs4,1)=mean([vSpotsXYZ{3}(i,1),vSpotsXYZ{4}(j,1)]);
                loc3Co4(n3Colocs4,2)=mean([vSpotsXYZ{3}(i,2),vSpotsXYZ{4}(j,2)]);
                loc3Co4(n3Colocs4,3)=mean([vSpotsXYZ{3}(i,3),vSpotsXYZ{4}(j,3)]);
                Stime34(n3Colocs4)=0;
                radii34(n3Colocs4)=SpotRadius;   
            else
                if ~isequal(Spot3Co4Logic(i),1)
                    Spot3Co4Logic(i)=0;
                end
              
                if ~isequal(Spot4Co3Logic(j),1)
                    Spot4Co3Logic(j)=0;
                end
               
            end    
        end
    end
end

%object based coloc-no double counting
if choice==2
    %%now build the triple & quadruple coloc by using an AND operation on the logical
    %%vectors for 1coloc2, 1coloc3 and 1coloc4, so we get a logical array 
    %%of spot 1 colocalized with 2, 3 and 4.
    
    Spot1co23Logic=Spot1Co2Logic & Spot1Co3Logic;
    Spot1co24Logic=Spot1Co2Logic & Spot1Co4Logic;
    Spot1co34Logic=Spot1Co3Logic & Spot1Co4Logic;
    
    Spot2co13Logic=Spot2Co1Logic & Spot2Co3Logic;
    Spot2co14Logic=Spot2Co1Logic & Spot2Co4Logic;
    Spot2co34Logic=Spot2Co3Logic & Spot2Co4Logic;
    
    Spot3co12Logic=Spot3Co1Logic & Spot3Co2Logic;
    Spot3co14Logic=Spot3Co1Logic & Spot3Co4Logic;
    Spot3co24Logic=Spot3Co2Logic & Spot3Co4Logic;

    Spot4co12Logic=Spot4Co1Logic & Spot4Co2Logic;
    Spot4co13Logic=Spot4Co1Logic & Spot4Co3Logic;    
    Spot4co23Logic=Spot4Co2Logic & Spot4Co3Logic;    

    Spot1co234Logic=Spot1Co2Logic & Spot1Co3Logic & Spot1Co4Logic;        
    
    [nlog12,p]=size(find(Spot1Co2Logic==1));
    [nlog21,p]=size(find(Spot2Co1Logic==1));    
    [nlog13,p]=size(find(Spot1Co3Logic==1));
    [nlog31,p]=size(find(Spot3Co1Logic==1));
    [nlog14,p]=size(find(Spot1Co4Logic==1));
    [nlog41,p]=size(find(Spot4Co1Logic==1));
    [nlog23,p]=size(find(Spot2Co3Logic==1));
    [nlog32,p]=size(find(Spot3Co2Logic==1));
    [nlog24,p]=size(find(Spot2Co4Logic==1));
    [nlog42,p]=size(find(Spot4Co2Logic==1));
    [nlog34,p]=size(find(Spot3Co4Logic==1));
    [nlog43,p]=size(find(Spot4Co3Logic==1));
    
    [nlog123,p]=size(find(Spot1co23Logic==1));
    [nlog124,p]=size(find(Spot1co24Logic==1));
    [nlog134,p]=size(find(Spot1co34Logic==1));

    [nlog213,p]=size(find(Spot2co13Logic==1));
    [nlog214,p]=size(find(Spot2co14Logic==1));
    [nlog234,p]=size(find(Spot2co34Logic==1));
    
    [nlog312,p]=size(find(Spot3co12Logic==1));
    [nlog314,p]=size(find(Spot3co14Logic==1));
    [nlog324,p]=size(find(Spot3co24Logic==1));

    [nlog1234,p]=size(find(Spot1co234Logic==1));
    

    %use logical vectors for getting location arrays of colocalized spots
    XYZ1w2=vSpotsXYZ{1}(Spot1Co2Logic==1,:);
    XYZ2w1=vSpotsXYZ{2}(Spot2Co1Logic==1,:);
    XYZ1w3=vSpotsXYZ{1}(Spot1Co3Logic==1,:);
    XYZ3w1=vSpotsXYZ{3}(Spot3Co1Logic==1,:);
    XYZ1w4=vSpotsXYZ{1}(Spot1Co4Logic==1,:);
    XYZ4w1=vSpotsXYZ{4}(Spot4Co1Logic==1,:);
    XYZ2w3=vSpotsXYZ{2}(Spot2Co3Logic==1,:);
    XYZ3w2=vSpotsXYZ{3}(Spot3Co2Logic==1,:);
    XYZ2w4=vSpotsXYZ{2}(Spot2Co4Logic==1,:);
    XYZ4w2=vSpotsXYZ{4}(Spot4Co2Logic==1,:);
    XYZ3w4=vSpotsXYZ{3}(Spot3Co4Logic==1,:);
    XYZ4w3=vSpotsXYZ{4}(Spot4Co3Logic==1,:);
    
    XYZ1co23=vSpotsXYZ{1}(Spot1co23Logic==1,:);
    XYZ1co24=vSpotsXYZ{1}(Spot1co24Logic==1,:);
    XYZ1co34=vSpotsXYZ{1}(Spot1co34Logic==1,:);

    XYZ2co13=vSpotsXYZ{2}(Spot2co13Logic==1,:);
    XYZ2co14=vSpotsXYZ{2}(Spot2co14Logic==1,:);
    XYZ2co34=vSpotsXYZ{2}(Spot2co34Logic==1,:);

    XYZ3co12=vSpotsXYZ{3}(Spot3co12Logic==1,:);
    XYZ3co14=vSpotsXYZ{3}(Spot3co14Logic==1,:);
    XYZ3co24=vSpotsXYZ{3}(Spot3co24Logic==1,:);

    XYZ1co234=vSpotsXYZ{1}(Spot1co234Logic==1,:);

    %get size of location arrays for building time and radii arrays
    [n1w2,p]=size(XYZ1w2);
    [n2w1,p]=size(XYZ2w1);
    [n1w3,p]=size(XYZ1w3);
    [n3w1,p]=size(XYZ3w1);
    [n1w4,p]=size(XYZ1w4);
    [n4w1,p]=size(XYZ4w1);
    [n2w3,p]=size(XYZ2w3);
    [n3w2,p]=size(XYZ3w2);
    [n2w4,p]=size(XYZ2w4);
    [n4w2,p]=size(XYZ4w2);
    [n3w4,p]=size(XYZ3w4);
    [n4w3,p]=size(XYZ4w3);
    
    [n123,p]=size(XYZ1co23);
    [n124,p]=size(XYZ1co24);
    [n134,p]=size(XYZ1co34);
    [n213,p]=size(XYZ2co13);
    [n214,p]=size(XYZ2co14);
    [n234,p]=size(XYZ2co34);
    [n312,p]=size(XYZ3co12);
    [n314,p]=size(XYZ3co14);
    [n324,p]=size(XYZ3co24);

    [n1234,p]=size(XYZ1co234);

    %build time arrays
    t1w2=zeros(n1w2,1);
    t2w1=zeros(n2w1,1);
    t1w3=zeros(n1w3,1);
    t3w1=zeros(n3w1,1);
    t1w4=zeros(n1w4,1);
    t4w1=zeros(n4w1,1);
    t2w3=zeros(n2w3,1);
    t3w2=zeros(n3w2,1);
    t2w4=zeros(n2w4,1);
    t4w2=zeros(n4w2,1);
    t3w4=zeros(n3w4,1);
    t4w3=zeros(n4w3,1);
    t123=zeros(n123,1);
    t124=zeros(n124,1);
    t134=zeros(n134,1);
    t213=zeros(n213,1);
    t214=zeros(n214,1);
    t234=zeros(n234,1);
    t312=zeros(n312,1);
    t314=zeros(n314,1);
    t324=zeros(n324,1);
    t1234=zeros(n1234,1);

    %build radius arrays
    r1w2=ones(n1w2,1)*vSpotsRadius{1}(1,1);
    r2w1=ones(n2w1,1)*vSpotsRadius{2}(1,1);
    r1w3=ones(n1w3,1)*vSpotsRadius{1}(1,1);
    r3w1=ones(n3w1,1)*vSpotsRadius{3}(1,1);
    r1w4=ones(n1w4,1)*vSpotsRadius{1}(1,1);
    r4w1=ones(n4w1,1)*vSpotsRadius{4}(1,1);
    r2w3=ones(n2w3,1)*vSpotsRadius{2}(1,1);
    r3w2=ones(n3w2,1)*vSpotsRadius{3}(1,1);
    r2w4=ones(n2w4,1)*vSpotsRadius{2}(1,1);
    r4w2=ones(n4w2,1)*vSpotsRadius{4}(1,1);
    r3w4=ones(n3w4,1)*vSpotsRadius{3}(1,1);
    r4w3=ones(n4w3,1)*vSpotsRadius{4}(1,1); 
    r123=ones(n123,1)*vSpotsRadius{1}(1,1);
    r124=ones(n124,1)*vSpotsRadius{1}(1,1);
    r134=ones(n134,1)*vSpotsRadius{1}(1,1);
    r213=ones(n213,1)*vSpotsRadius{2}(1,1);
    r214=ones(n214,1)*vSpotsRadius{2}(1,1);
    r234=ones(n234,1)*vSpotsRadius{2}(1,1);
    r312=ones(n312,1)*vSpotsRadius{3}(1,1);
    r314=ones(n314,1)*vSpotsRadius{3}(1,1);
    r324=ones(n324,1)*vSpotsRadius{3}(1,1);
    r1234=ones(n1234,1)*vSpotsRadius{1}(1,1);

%% create and add spots for user verification of coloc results
    if ~isempty(XYZ1w2)
        vNewSpots1 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots1.Set(XYZ1w2, t1w2, r1w2);
        vNewSpots1.SetName([spotName{1},' colocated with ',spotName{2}]);
        %set color to yellow
        vRGBA = 65279;
        vNewSpots1.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots1, -1);
    end

    if ~isempty(XYZ2w1)
        vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots2.Set(XYZ2w1, t2w1, r2w1);
        vNewSpots2.SetName([spotName{2},' colocated with ',spotName{1}]);
        %set color to light blue
        vRGBA = 16776704;
        vNewSpots2.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots2, -1);
    end

    if ~isempty(XYZ1w3)
        vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots2.Set(XYZ1w3, t1w3, r1w3);
        vNewSpots2.SetName([spotName{1},' colocated with ',spotName{3}]);
        %set color to light blue
        vRGBA = 16776704;
        vNewSpots2.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots2, -1);
    end

    if ~isempty(XYZ3w1)
        vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots2.Set(XYZ3w1, t3w1, r3w1);
        vNewSpots2.SetName([spotName{3},' colocated with ',spotName{1}]);
        %set color to purple
        vRGBA = 16711934;
        vNewSpots2.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots2, -1);
    end

    if ~isempty(XYZ1w4)
        vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots2.Set(XYZ1w4, t1w4, r1w4);
        vNewSpots2.SetName([spotName{1},' colocated with ',spotName{4}]);
        %set color to light blue
        vRGBA = 16776704;
        vNewSpots2.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots2, -1);
    end

    if ~isempty(XYZ4w1)
        vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots2.Set(XYZ4w1, t4w1, r4w1);
        vNewSpots2.SetName([spotName{4},' colocated with ',spotName{1}]);
        %set color to purple
        vRGBA = 16711934;
        vNewSpots2.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots2, -1);
    end
    
    if ~isempty(XYZ2w3)
        vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots2.Set(XYZ2w3, t2w3, r2w3);
        vNewSpots2.SetName([spotName{2},' colocated with ',spotName{3}]);
        %set color to light blue
        vRGBA = 16776704;
        vNewSpots2.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots2, -1);
    end
    
    if ~isempty(XYZ3w2)
        vNewSpots1 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots1.Set(XYZ3w2, t3w2, r3w2);
        vNewSpots1.SetName([spotName{3},' colocated with ',spotName{2}]);
        %set color to yellow
        vRGBA = 65279;
        vNewSpots1.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots1, -1);
    end

    if ~isempty(XYZ2w4)
        vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots2.Set(XYZ2w4, t2w4, r2w4);
        vNewSpots2.SetName([spotName{2},' colocated with ',spotName{4}]);
        %set color to purple
        vRGBA = 16711934;
        vNewSpots2.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots2, -1);
    end
    
    if ~isempty(XYZ4w2)
        vNewSpots1 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots1.Set(XYZ4w2, t4w2, r4w2);
        vNewSpots1.SetName([spotName{4},' colocated with ',spotName{2}]);
        %set color to yellow
        vRGBA = 65279;
        vNewSpots1.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots1, -1);
    end

    if ~isempty(XYZ3w4)
        vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots2.Set(XYZ3w4, t3w4, r3w4);
        vNewSpots2.SetName([spotName{3},' colocated with ',spotName{4}]);
        %set color to purple
        vRGBA = 16711934;
        vNewSpots2.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots2, -1);
    end
    
    if ~isempty(XYZ4w3)
        vNewSpots1 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots1.Set(XYZ4w3, t4w3, r4w3);
        vNewSpots1.SetName([spotName{4},' colocated with ',spotName{3}]);
        %set color to yellow
        vRGBA = 65279;
        vNewSpots1.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots1, -1);
    end
    
    if ~isempty(XYZ1co23)
        vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots2.Set(XYZ1co23, t123, r123);
        vNewSpots2.SetName([spotName{1},' colocated with ',spotName{2},' and ',spotName{3}]);
        %set color to WHITE
        vRGBA = 16776958;
        vNewSpots2.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots2, -1);
    end
    
    if ~isempty(XYZ1co24)
        vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots2.Set(XYZ1co24, t124, r124);
        vNewSpots2.SetName([spotName{1},' colocated with ',spotName{2},' and ',spotName{4}]);
        %set color to yellow
        vRGBA = 65279;
        vNewSpots2.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots2, -1);
    end

    if ~isempty(XYZ1co34)
        vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots2.Set(XYZ1co34, t134, r134);
        vNewSpots2.SetName([spotName{1},' colocated with ',spotName{3},' and ',spotName{4}]);
        %set color to WHITE
        vRGBA = 16776958;
        vNewSpots2.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots2, -1);
    end
    
    if ~isempty(XYZ2co13)
        vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots2.Set(XYZ2co13, t213, r213);
        vNewSpots2.SetName([spotName{2},' colocated with ',spotName{1}, ' and ',spotName{3}]);
        %set color to yellow
        vRGBA = 65279;
        vNewSpots2.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots2, -1);
    end
    
    if ~isempty(XYZ2co14)
        vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots2.Set(XYZ2co14, t214, r214);
        vNewSpots2.SetName([spotName{2},' colocated with ',spotName{1}, ' and ',spotName{4}]);
        %set color to WHITE
        vRGBA = 16776958;
        vNewSpots2.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots2, -1);
    end
    
    if ~isempty(XYZ2co34)
        vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots2.Set(XYZ2co34, t234, r234);
        vNewSpots2.SetName([spotName{2},' colocated with ',spotName{3}, ' and ',spotName{4}]);
        %set color to purple
        vRGBA = 16711934;
        vNewSpots2.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots2, -1);
    end
    
    if ~isempty(XYZ3co12)
        vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots2.Set(XYZ3co12, t312, r312);
        vNewSpots2.SetName([spotName{3},' colocated with ',spotName{1}, ' and ',spotName{2}]);
        %set color to WHITE
        vRGBA = 16776958;
        vNewSpots2.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots2, -1);
    end
    
    if ~isempty(XYZ3co14)
        vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots2.Set(XYZ3co14, t314, r314);
        vNewSpots2.SetName([spotName{3},' colocated with ',spotName{1}, ' and ',spotName{4}]);
        %set color to WHITE
        vRGBA = 16776958;
        vNewSpots2.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots2, -1);
    end    

    if ~isempty(XYZ3co24)
        vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots2.Set(XYZ3co24, t324, r324);
        vNewSpots2.SetName([spotName{3},' colocated with ',spotName{2}, ' and ',spotName{4}]);
        %set color to purple
        vRGBA = 16711934;
        vNewSpots2.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots2, -1);
    end
    
    if ~isempty(XYZ1co234)
        vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
        vNewSpots2.Set(XYZ1co234, t1234, r1234);
        vNewSpots2.SetName([spotName{1},'quadruple colocated']);
        %set color to purple
        vRGBA = 16711934;
        vNewSpots2.SetColorRGBA(vRGBA);
        vSurpassScene.AddChild(vNewSpots2, -1);
    end
    
%write data to file
    filename=char(vImarisApplication.GetCurrentFileName);
    filename2=regexprep(filename,'.ims','_coloc_spots_data.csv');
    filename2=regexprep(filename2,'.czi','_coloc_spots_data.csv');
    filename2=regexprep(filename2,'.lsm','_coloc_spots_data.csv');
    filename2=regexprep(filename2,'.tif','_coloc_spots_data.csv');
    filename2=regexprep(filename2,'.tiff','_coloc_spots_data.csv');
    [name3 pathname3]=uiputfile(filename2,'Choose where to Save CSV files');
    filename2 = [pathname3 name3];
    comma=',';
    fid = fopen(filename2, 'w');
    fprintf(fid, '%s', filename);
    fprintf(fid, '%s', comma);

    fprintf(fid, '%g', Nspots(1));
    fprintf(fid, '%s', comma);

    fprintf(fid, '%g', Nspots(2));
    fprintf(fid, '%s', comma);

    fprintf(fid, '%g', Nspots(3));
    fprintf(fid, '%s', comma);

    fprintf(fid, '%g', Nspots(4));
    fprintf(fid, '%s', comma);
    
    fprintf(fid, '%g', n1w2);
    fprintf(fid, '%s', comma);

    fprintf(fid, '%g', n2w1);
    fprintf(fid, '%s', comma);
    
    fprintf(fid, '%g', n1w3);
    fprintf(fid, '%s', comma);  

    fprintf(fid, '%g', n3w1);
    fprintf(fid, '%s', comma);

    fprintf(fid, '%g', n1w4);
    fprintf(fid, '%s', comma);  

    fprintf(fid, '%g', n4w1);
    fprintf(fid, '%s', comma);
    
    fprintf(fid, '%g', n2w3);
    fprintf(fid, '%s', comma);  

    fprintf(fid, '%g', n3w2);
    fprintf(fid, '%s', comma);

    fprintf(fid, '%g', n2w4);
    fprintf(fid, '%s', comma);  

    fprintf(fid, '%g', n4w2);
    fprintf(fid, '%s', comma);    

    fprintf(fid, '%g', n3w4);
    fprintf(fid, '%s', comma);  

    fprintf(fid, '%g', n4w3);
    fprintf(fid, '%s', comma);
    
    fprintf(fid, '%g', n123);
    fprintf(fid, '%s', comma);  

    fprintf(fid, '%g', n124);
    fprintf(fid, '%s', comma);  
    
    fprintf(fid, '%g', n134);
    fprintf(fid, '%s', comma);  
    
    fprintf(fid, '%g', n213);
    fprintf(fid, '%s', comma);
    
    fprintf(fid, '%g', n214);
    fprintf(fid, '%s', comma);
    
    fprintf(fid, '%g', n234);
    fprintf(fid, '%s', comma);
    
    fprintf(fid, '%g', n312);
    fprintf(fid, '%s', comma);  

    fprintf(fid, '%g', n314);
    fprintf(fid, '%s', comma); 
    
    fprintf(fid, '%g', n324);
    fprintf(fid, '%s', comma); 

    fprintf(fid, '%g', n1234);
    fprintf(fid, '%s', comma);
    
    fprintf(fid, '%s', [spotName{1},'-',spotName{2},'-',spotName{3},'-',spotName{4}]);
    fprintf(fid, '%s', comma);

    fprintf(fid, '%s\n', ' ');


    fclose(fid);
    
    
    for i=1:nSpotsObj
        newFileName{i}=['_' spotName{i} '.csv'];
    end    
    filename3=regexprep(filename2,'_coloc_spots_data.csv',newFileName{1});
    filename4=regexprep(filename2,'_coloc_spots_data.csv',newFileName{2});
    filename5=regexprep(filename2,'_coloc_spots_data.csv',newFileName{3});
    filename6=regexprep(filename2,'_coloc_spots_data.csv',newFileName{4});

    if nSpotsObj>4
        filename7=regexprep(filename2,'_coloc_spots_data.csv',newFileName{5});
    end
    if nSpotsObj>5
        filename8=regexprep(filename2,'_coloc_spots_data.csv',newFileName{6});
    end
    if nSpotsObj>6
        filename9=regexprep(filename2,'_coloc_spots_data.csv',newFileName{7});
    end   
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    xOffset=0;
    yOffset=0;
    lScale=1;
    bScale=1;

    if ~isempty(vSpotsXYZ{4})
        vSpotsXYZ{4}(:,1)=vSpotsXYZ{4}(:,1)-xOffset; 
        vSpotsXYZ{4}(:,2)=vSpotsXYZ{4}(:,2)-yOffset;
        vSpotsXYZ{4}(:,1)=vSpotsXYZ{4}(:,1)*lScale; 
        vSpotsXYZ{4}(:,2)=vSpotsXYZ{4}(:,2)*bScale;
    end
    
    if ~isempty(vSpotsXYZ{3})
        vSpotsXYZ{3}(:,1)=vSpotsXYZ{3}(:,1)-xOffset; 
        vSpotsXYZ{3}(:,2)=vSpotsXYZ{3}(:,2)-yOffset;
        vSpotsXYZ{3}(:,1)=vSpotsXYZ{3}(:,1)*lScale; 
        vSpotsXYZ{3}(:,2)=vSpotsXYZ{3}(:,2)*bScale;
    end
    if ~isempty(vSpotsXYZ{2})
        vSpotsXYZ{2}(:,1)=vSpotsXYZ{2}(:,1)-xOffset; 
        vSpotsXYZ{2}(:,2)=vSpotsXYZ{2}(:,2)-yOffset;
        vSpotsXYZ{2}(:,1)=vSpotsXYZ{2}(:,1)*lScale; 
        vSpotsXYZ{2}(:,2)=vSpotsXYZ{2}(:,2)*bScale;
    end
    if ~isempty(vSpotsXYZ{1})
        vSpotsXYZ{1}(:,1)=vSpotsXYZ{1}(:,1)-xOffset; 
        vSpotsXYZ{1}(:,2)=vSpotsXYZ{1}(:,2)-yOffset;
        vSpotsXYZ{1}(:,1)=vSpotsXYZ{1}(:,1)*lScale; 
        vSpotsXYZ{1}(:,2)=vSpotsXYZ{1}(:,2)*bScale;
    end
    if ~isempty(XYZ1w2)
        XYZ1w2(:,1)=XYZ1w2(:,1)-xOffset; 
        XYZ1w2(:,2)=XYZ1w2(:,2)-yOffset;
        XYZ1w2(:,1)=XYZ1w2(:,1)*lScale; 
        XYZ1w2(:,2)=XYZ1w2(:,2)*bScale;
    end
    if ~isempty(XYZ2w1)
        XYZ2w1(:,1)=XYZ2w1(:,1)-xOffset; 
        XYZ2w1(:,2)=XYZ2w1(:,2)-yOffset;
        XYZ2w1(:,1)=XYZ2w1(:,1)*lScale; 
        XYZ2w1(:,2)=XYZ2w1(:,2)*bScale;
    end
    if ~isempty(XYZ1w3)
        XYZ1w3(:,1)=XYZ1w3(:,1)-xOffset; 
        XYZ1w3(:,2)=XYZ1w3(:,2)-yOffset;
        XYZ1w3(:,1)=XYZ1w3(:,1)*lScale; 
        XYZ1w3(:,2)=XYZ1w3(:,2)*bScale;
    end
    if ~isempty(XYZ3w1)
        XYZ3w1(:,1)=XYZ3w1(:,1)-xOffset; 
        XYZ3w1(:,2)=XYZ3w1(:,2)-yOffset;
        XYZ3w1(:,1)=XYZ3w1(:,1)*lScale; 
        XYZ3w1(:,2)=XYZ3w1(:,2)*bScale;
    end
    if ~isempty(XYZ1w4)
        XYZ1w4(:,1)=XYZ1w4(:,1)-xOffset; 
        XYZ1w4(:,2)=XYZ1w4(:,2)-yOffset;
        XYZ1w4(:,1)=XYZ1w4(:,1)*lScale; 
        XYZ1w4(:,2)=XYZ1w4(:,2)*bScale;
    end
    if ~isempty(XYZ4w1)
        XYZ4w1(:,1)=XYZ4w1(:,1)-xOffset; 
        XYZ4w1(:,2)=XYZ4w1(:,2)-yOffset;
        XYZ4w1(:,1)=XYZ4w1(:,1)*lScale; 
        XYZ4w1(:,2)=XYZ4w1(:,2)*bScale;
    end
    if ~isempty(XYZ2w3)
        XYZ2w3(:,1)=XYZ2w3(:,1)-xOffset; 
        XYZ2w3(:,2)=XYZ2w3(:,2)-yOffset;
        XYZ2w3(:,1)=XYZ2w3(:,1)*lScale; 
        XYZ2w3(:,2)=XYZ2w3(:,2)*bScale;
    end
    if ~isempty(XYZ3w2)
        XYZ3w2(:,1)=XYZ3w2(:,1)-xOffset; 
        XYZ3w2(:,2)=XYZ3w2(:,2)-yOffset;
        XYZ3w2(:,1)=XYZ3w2(:,1)*lScale; 
        XYZ3w2(:,2)=XYZ3w2(:,2)*bScale;
    end
    if ~isempty(XYZ2w4)
        XYZ2w4(:,1)=XYZ2w4(:,1)-xOffset; 
        XYZ2w4(:,2)=XYZ2w4(:,2)-yOffset;
        XYZ2w4(:,1)=XYZ2w4(:,1)*lScale; 
        XYZ2w4(:,2)=XYZ2w4(:,2)*bScale;
    end
    if ~isempty(XYZ4w2)
        XYZ4w2(:,1)=XYZ4w2(:,1)-xOffset; 
        XYZ4w2(:,2)=XYZ4w2(:,2)-yOffset;
        XYZ4w2(:,1)=XYZ4w2(:,1)*lScale; 
        XYZ4w2(:,2)=XYZ4w2(:,2)*bScale;
    end
    if ~isempty(XYZ3w4)
        XYZ3w4(:,1)=XYZ3w4(:,1)-xOffset; 
        XYZ3w4(:,2)=XYZ3w4(:,2)-yOffset;
        XYZ3w4(:,1)=XYZ3w4(:,1)*lScale; 
        XYZ3w4(:,2)=XYZ3w4(:,2)*bScale;
    end
    if ~isempty(XYZ4w3)
        XYZ4w3(:,1)=XYZ4w3(:,1)-xOffset; 
        XYZ4w3(:,2)=XYZ4w3(:,2)-yOffset;
        XYZ4w3(:,1)=XYZ4w3(:,1)*lScale; 
        XYZ4w3(:,2)=XYZ4w3(:,2)*bScale;
    end
    if ~isempty(XYZ1co23)
        XYZ1co23(:,1)=XYZ1co23(:,1)-xOffset; 
        XYZ1co23(:,2)=XYZ1co23(:,2)-yOffset;
        XYZ1co23(:,1)=XYZ1co23(:,1)*lScale; 
        XYZ1co23(:,2)=XYZ1co23(:,2)*bScale;
    end
    if ~isempty(XYZ1co24)
        XYZ1co24(:,1)=XYZ1co24(:,1)-xOffset; 
        XYZ1co24(:,2)=XYZ1co24(:,2)-yOffset;
        XYZ1co24(:,1)=XYZ1co24(:,1)*lScale; 
        XYZ1co24(:,2)=XYZ1co24(:,2)*bScale;
    end
    if ~isempty(XYZ1co34)
        XYZ1co34(:,1)=XYZ1co34(:,1)-xOffset; 
        XYZ1co34(:,2)=XYZ1co34(:,2)-yOffset;
        XYZ1co34(:,1)=XYZ1co34(:,1)*lScale; 
        XYZ1co34(:,2)=XYZ1co34(:,2)*bScale;
    end
    if ~isempty(XYZ2co13)
        XYZ2co13(:,1)=XYZ2co13(:,1)-xOffset; 
        XYZ2co13(:,2)=XYZ2co13(:,2)-yOffset;
        XYZ2co13(:,1)=XYZ2co13(:,1)*lScale; 
        XYZ2co13(:,2)=XYZ2co13(:,2)*bScale;
    end
    if ~isempty(XYZ2co14)
        XYZ2co14(:,1)=XYZ2co14(:,1)-xOffset; 
        XYZ2co14(:,2)=XYZ2co14(:,2)-yOffset;
        XYZ2co14(:,1)=XYZ2co14(:,1)*lScale; 
        XYZ2co14(:,2)=XYZ2co14(:,2)*bScale;
    end
    if ~isempty(XYZ2co34)
        XYZ2co34(:,1)=XYZ2co34(:,1)-xOffset; 
        XYZ2co34(:,2)=XYZ2co34(:,2)-yOffset;
        XYZ2co34(:,1)=XYZ2co34(:,1)*lScale; 
        XYZ2co34(:,2)=XYZ2co34(:,2)*bScale;
    end
    if ~isempty(XYZ3co12)
        XYZ3co12(:,1)=XYZ3co12(:,1)-xOffset; 
        XYZ3co12(:,2)=XYZ3co12(:,2)-yOffset;
        XYZ3co12(:,1)=XYZ3co12(:,1)*lScale; 
        XYZ3co12(:,2)=XYZ3co12(:,2)*bScale;
    end
    if ~isempty(XYZ3co14)
        XYZ3co14(:,1)=XYZ3co14(:,1)-xOffset; 
        XYZ3co14(:,2)=XYZ3co14(:,2)-yOffset;
        XYZ3co14(:,1)=XYZ3co14(:,1)*lScale; 
        XYZ3co14(:,2)=XYZ3co14(:,2)*bScale;
    end    
    if ~isempty(XYZ3co24)
        XYZ3co24(:,1)=XYZ3co24(:,1)-xOffset; 
        XYZ3co24(:,2)=XYZ3co24(:,2)-yOffset;
        XYZ3co24(:,1)=XYZ3co24(:,1)*lScale; 
        XYZ3co24(:,2)=XYZ3co24(:,2)*bScale;
    end    
    if ~isempty(XYZ1co234)
        XYZ1co234(:,1)=XYZ1co234(:,1)-xOffset; 
        XYZ1co234(:,2)=XYZ1co234(:,2)-yOffset;
        XYZ1co234(:,1)=XYZ1co234(:,1)*lScale; 
        XYZ1co234(:,2)=XYZ1co234(:,2)*bScale;
    end    
    
    catMatrix=zeros([6 3]);
    catMatrix(2,1:2)=ccPoint(1:2);
    catMatrix(3,1:2)=bPoint(1:2);
    catMatrix(4,1:2)=lPoint(1:2);
    catMatrix(5,1:2)=tPoint(1:2);
    catMatrix(6,1:2)=[Rotation,0];
    catMatrix(:,1)=catMatrix(:,1)-xOffset;
    catMatrix(:,2)=catMatrix(:,2)-yOffset;
    catMatrix(:,1)=catMatrix(:,1)*lScale;
    catMatrix(:,2)=catMatrix(:,2)*bScale;
    catMatrix(1,:)=0;
    
    
    
    if ~isempty(vSpotsXYZ{4})
        vSpotsXYZ{4}=cat(1,vSpotsXYZ{4},catMatrix);
    else
        vSpotsXYZ{4}=catMatrix;
    end    
    if ~isempty(vSpotsXYZ{3})
        vSpotsXYZ{3}=cat(1,vSpotsXYZ{3},catMatrix);
    else
        vSpotsXYZ{3}=catMatrix;
    end
    if ~isempty(vSpotsXYZ{2})
        vSpotsXYZ{2}=cat(1,vSpotsXYZ{2},catMatrix);
    else
        vSpotsXYZ{2}=catMatrix;    
    end
    if ~isempty(vSpotsXYZ{1})
        vSpotsXYZ{1}=cat(1,vSpotsXYZ{1},catMatrix);
    else
        vSpotsXYZ{1}=catMatrix;    
    end
    if ~isempty(XYZ1w2)
        XYZ1w2=cat(1,XYZ1w2,catMatrix);
    else
        XYZ1w2=catMatrix;
    end
    if ~isempty(XYZ2w1)
        XYZ2w1=cat(1,XYZ2w1,catMatrix);
    else
        XYZ2w1=catMatrix;
    end
    if ~isempty(XYZ1w3)
        XYZ1w3=cat(1,XYZ1w3,catMatrix);
    else
        XYZ1w3=catMatrix;
    end
    if ~isempty(XYZ3w1)
        XYZ3w1=cat(1,XYZ3w1,catMatrix);
    else
        XYZ3w1=catMatrix;
    end
    if ~isempty(XYZ1w4)
        XYZ1w4=cat(1,XYZ1w4,catMatrix);
    else
        XYZ1w4=catMatrix;
    end
    if ~isempty(XYZ4w1)
        XYZ4w1=cat(1,XYZ4w1,catMatrix);
    else
        XYZ4w1=catMatrix;
    end
    if ~isempty(XYZ2w3)
        XYZ2w3=cat(1,XYZ2w3,catMatrix);
    else
        XYZ2w3=catMatrix;
    end
    if ~isempty(XYZ3w2)
        XYZ3w2=cat(1,XYZ3w2,catMatrix);
    else
        XYZ3w2=catMatrix;
    end
    if ~isempty(XYZ2w4)
        XYZ2w4=cat(1,XYZ2w4,catMatrix);
    else
        XYZ2w4=catMatrix;
    end
    if ~isempty(XYZ4w2)
        XYZ4w2=cat(1,XYZ4w2,catMatrix);
    else
        XYZ4w2=catMatrix;
    end
    if ~isempty(XYZ3w4)
        XYZ3w4=cat(1,XYZ3w4,catMatrix);
    else
        XYZ3w4=catMatrix;
    end
    if ~isempty(XYZ4w3)
        XYZ4w3=cat(1,XYZ4w3,catMatrix);
    else
        XYZ4w3=catMatrix;
    end
    if ~isempty(XYZ1co23)
        XYZ1co23=cat(1,XYZ1co23,catMatrix);
    else
        XYZ1co23=catMatrix;
    end
    if ~isempty(XYZ1co24)
        XYZ1co24=cat(1,XYZ1co24,catMatrix);
    else
        XYZ1co24=catMatrix;
    end
    if ~isempty(XYZ1co34)
        XYZ1co34=cat(1,XYZ1co34,catMatrix);
    else
        XYZ1co34=catMatrix;
    end
    if ~isempty(XYZ2co13)
        XYZ2co13=cat(1,XYZ2co13,catMatrix);
    else
        XYZ2co13=catMatrix;
    end
    if ~isempty(XYZ2co14)
        XYZ2co14=cat(1,XYZ2co14,catMatrix);
    else
        XYZ2co14=catMatrix;
    end
    if ~isempty(XYZ2co34)
        XYZ2co34=cat(1,XYZ2co34,catMatrix);
    else
        XYZ2co34=catMatrix;
    end
    if ~isempty(XYZ3co12)
        XYZ3co12=cat(1,XYZ3co12,catMatrix);
    else
        XYZ3co12=catMatrix;
    end
    if ~isempty(XYZ3co14)
        XYZ3co14=cat(1,XYZ3co14,catMatrix);
    else
        XYZ3co14=catMatrix;
    end
    if ~isempty(XYZ3co24)
        XYZ3co24=cat(1,XYZ3co24,catMatrix);
    else
        XYZ3co24=catMatrix;
    end    
    headersForCSVs='X coordinate (scaled to Reference SC), Y coordinate (scaled to Reference SC)';
    
    
    [len1,width1]=size(vSpotsXYZ{1});
    if len1>6 %indicating vSpotsXYZ{1} was nonempty before concatenation.
        fid = fopen(filename3, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len1
            fprintf(fid, '%g', vSpotsXYZ{1}(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', vSpotsXYZ{1}(i,2));
    %         fprintf(fid, '%s', comma);
    %         fprintf(fid, '%g', vSpotsXYZ{1}(i,3));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end    
    
    [len2,width2]=size(vSpotsXYZ{2});
    if len2>6
        fid = fopen(filename4, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len2
            fprintf(fid, '%g', vSpotsXYZ{2}(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', vSpotsXYZ{2}(i,2));
    %         fprintf(fid, '%s', comma);
    %         fprintf(fid, '%g', vSpotsXYZ{2}(i,3));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end    
    
    [len3,width3]=size(vSpotsXYZ{3});
    if len3>6
        fid = fopen(filename5, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len3
            fprintf(fid, '%g', vSpotsXYZ{3}(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', vSpotsXYZ{3}(i,2));
    %         fprintf(fid, '%s', comma);
    %         fprintf(fid, '%g', vSpotsXYZ{3}(i,3));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end

    [len4,width4]=size(vSpotsXYZ{4});
    if len4>6
        fid = fopen(filename6, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len4
            fprintf(fid, '%g', vSpotsXYZ{4}(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', vSpotsXYZ{4}(i,2));
    %         fprintf(fid, '%s', comma);
    %         fprintf(fid, '%g', vSpotsXYZ{3}(i,3));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end
    
    [len5,width5]=size(XYZ1w2);
    if len5>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
            spotName{1} '_w_' spotName{2} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len5
            fprintf(fid, '%g', XYZ1w2(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ1w2(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end
    
    [len6,width6]=size(XYZ2w1);
    if len6>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
            spotName{2} '_w_' spotName{1} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len6
            fprintf(fid, '%g', XYZ2w1(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ2w1(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end

    [len7,width7]=size(XYZ1w3);
    if len7>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
            spotName{1} '_w_' spotName{3} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len7
            fprintf(fid, '%g', XYZ1w3(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ1w3(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end
    
    [len8,width8]=size(XYZ3w1);
    if len8>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
            spotName{3} '_w_' spotName{1} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len8
            fprintf(fid, '%g', XYZ3w1(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ3w1(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end

    [len9,width9]=size(XYZ1w4);
    if len9>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
            spotName{1} '_w_' spotName{4} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len9
            fprintf(fid, '%g', XYZ1w4(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ1w4(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end
    
    [len10,width10]=size(XYZ4w1);
    if len10>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
            spotName{4} '_w_' spotName{1} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len10
            fprintf(fid, '%g', XYZ4w1(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ4w1(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end
    
    [len11,width11]=size(XYZ2w3);
    if len11>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
            spotName{2} '_w_' spotName{3} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len11
            fprintf(fid, '%g', XYZ2w3(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ2w3(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end
    
    [len12,width12]=size(XYZ3w2);
    if len12>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
            spotName{3} '_w_' spotName{2} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len12
            fprintf(fid, '%g', XYZ3w2(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ3w2(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end
 
    [len13,width13]=size(XYZ2w4);
    if len13>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
            spotName{2} '_w_' spotName{4} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len13
            fprintf(fid, '%g', XYZ2w4(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ2w4(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end
    
    [len14,width14]=size(XYZ4w2);
    if len14>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
            spotName{4} '_w_' spotName{2} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len14
            fprintf(fid, '%g', XYZ4w2(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ4w2(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end
    
    [len15,width15]=size(XYZ3w4);
    if len15>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
            spotName{3} '_w_' spotName{4} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len15
            fprintf(fid, '%g', XYZ3w4(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ3w4(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end
    
    [len16,width16]=size(XYZ4w3);
    if len16>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
            spotName{4} '_w_' spotName{3} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len16
            fprintf(fid, '%g', XYZ4w3(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ4w3(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end
    
    [len17,width17]=size(XYZ1co23);
    if len17>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
            spotName{1} '_w_' spotName{2} '&' spotName{3} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len17
            fprintf(fid, '%g', XYZ1co23(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ1co23(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end

    [len18,width18]=size(XYZ1co24);
    if len18>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
                 spotName{1} '_w_' spotName{2} '&' spotName{4} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len18
            fprintf(fid, '%g', XYZ1co24(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ1co24(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end
 
    [len19,width19]=size(XYZ1co34);
    if len19>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
                 spotName{1} '_w_' spotName{3} '&' spotName{4} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len19
            fprintf(fid, '%g', XYZ1co34(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ1co34(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end
    
    [len20,width20]=size(XYZ2co13);
    if len20>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
            spotName{2} '_w_' spotName{1} '&' spotName{3} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len20
            fprintf(fid, '%g', XYZ2co13(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ2co13(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end

    [len21,width21]=size(XYZ2co14);
    if len21>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
            spotName{2} '_w_' spotName{1} '&' spotName{4} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len21
            fprintf(fid, '%g', XYZ2co14(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ2co14(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end

    [len22,width22]=size(XYZ2co34);
    if len22>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
            spotName{2} '_w_' spotName{3} '&' spotName{4} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len22
            fprintf(fid, '%g', XYZ2co34(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ2co34(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end
    
    [len23,width23]=size(XYZ3co12);
    if len23>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
            spotName{3} '_w_' spotName{1} '&' spotName{2} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len23
            fprintf(fid, '%g', XYZ3co12(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ3co12(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end
    
    [len24,width24]=size(XYZ3co14);
    if len24>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
            spotName{3} '_w_' spotName{1} '&' spotName{4} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len24
            fprintf(fid, '%g', XYZ3co14(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ3co14(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end    

    [len25,width25]=size(XYZ3co24);
    if len25>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
            spotName{3} '_w_' spotName{2} '&' spotName{4} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len25
            fprintf(fid, '%g', XYZ3co24(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ3co24(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end  

    [len26,width26]=size(XYZ1co234);
    if len26>6
        thisName=regexprep(filename2,'_coloc_spots_data.csv',['_' ...
            spotName{1} '_w_' spotName{2} '&' spotName{3} '&' spotName{4} '.csv']);
        fid = fopen(thisName, 'w');
        fprintf(fid, '%s', headersForCSVs);
        fprintf(fid, '%s\n', ' ');
        for i=1:len26
            fprintf(fid, '%g', XYZ1co234(i,1));
            fprintf(fid, '%s', comma);
            fprintf(fid, '%g', XYZ1co234(i,2));
            fprintf(fid, '%s\n', ' ');
        end    
        fclose(fid);
    end
    
    stopHere2Debug=3;
    
end    
%% position based coloc---risk of double counting events
% if choice==1
%     %Find Coloc for Spot12 and Spot23
%     n12Colocs23=0;
%     n123Colocs=0;
%     loc12Co23=[];
%     if ~isempty(vSpotsXYZ{1}) || ~isempty(vSpotsXYZ{3})
%         [Nspots12,dim]=size(loc1Co2);
%         [Nspots23,dim]=size(loc2Co3);
%         [Nspots13,dim]=size(loc1Co3);
%         for i=1:Nspots12
%             for j=1:Nspots23
%                 Dx=loc1Co2(i,1)-loc2Co3(j,1);
%                 Dy=loc1Co2(i,2)-loc2Co3(j,2);
%                 Dz=loc1Co2(i,3)-loc2Co3(j,3);
%                 Dist12=(Dx^2)+(Dy^2)+(Dz^2);
%                 if Dist12<=MaxDist
%                     n12Colocs23=n12Colocs23+1;
%                     pos12Co23(n12Colocs23,1)=mean([loc1Co2(i,1),loc2Co3(j,1)]);
%                     pos12Co23(n12Colocs23,2)=mean([loc1Co2(i,2),loc2Co3(j,2)]);
%                     pos12Co23(n12Colocs23,3)=mean([loc1Co2(i,3),loc2Co3(j,3)]);
%                     for k=1:Nspots13
%                         D2x=pos12Co23(n12Colocs23,1)-loc1Co3(k,1);
%                         D2y=pos12Co23(n12Colocs23,2)-loc1Co3(k,2);
%                         D2z=pos12Co23(n12Colocs23,3)-loc1Co3(k,3);
%                         Dist123=(D2x^2)+(D2y^2)+(D2z^2);
%                         if Dist123<=MaxDist
%                             n123Colocs=n123Colocs+1;
%                             loc12Co23(n123Colocs,1)=mean([pos12Co23(n12Colocs23,1),loc1Co3(k,1)]);
%                             loc12Co23(n123Colocs,2)=mean([pos12Co23(n12Colocs23,2),loc1Co3(k,2)]);
%                             loc12Co23(n123Colocs,3)=mean([pos12Co23(n12Colocs23,3),loc1Co3(k,3)]);
%                             Stime1223(n123Colocs)=0;
%                             radii1223(n123Colocs)=SpotRadius;
%                         end
%                     end    
%                 end    
%             end
%         end
%     end  
%     %% create and add spots for user verification of coloc results
%     if ~isempty(loc1Co2)
%         vNewSpots1 = vImarisApplication.GetFactory.CreateSpots;
%         vNewSpots1.Set(loc1Co2, Stime12, radii12);
%         vNewSpots1.SetName([spotName{1},'-',spotName{2},' colocated']);
%         %set color to yellow
%         vRGBA = 65279;
%         vNewSpots1.SetColorRGBA(vRGBA);
%         vSurpassScene.AddChild(vNewSpots1, -1);
%     end
% 
%     if ~isempty(loc2Co3)
%         vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
%         vNewSpots2.Set(loc2Co3, Stime23, radii23);
%         vNewSpots2.SetName([spotName{2},'-',spotName{3},' colocated']);
%         %set color to light blue
%         vRGBA = 16776704;
%         vNewSpots2.SetColorRGBA(vRGBA);
%         vSurpassScene.AddChild(vNewSpots2, -1);
%     end
% 
%     if ~isempty(loc1Co3)
%         vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
%         vNewSpots2.Set(loc1Co3, Stime13, radii13);
%         vNewSpots2.SetName([spotName{1},'-',spotName{3},' colocated']);
%         %set color to purple
%         vRGBA = 16711934;
%         vNewSpots2.SetColorRGBA(vRGBA);
%         vSurpassScene.AddChild(vNewSpots2, -1);
%     end
% 
%     if ~isempty(loc12Co23)
%         vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
%         vNewSpots2.Set(loc12Co23, Stime1223, radii1223);
%         vNewSpots2.SetName([spotName{1},'-',spotName{2},'-',spotName{3},' colocated']);
%         %set color to WHITE
%         vRGBA = 16776958;
%         vNewSpots2.SetColorRGBA(vRGBA);
%         vSurpassScene.AddChild(vNewSpots2, -1);
%     end
% 
% %write data to file
%     filename=char(vImarisApplication.GetCurrentFileName);
%     filename2=regexprep(filename,'.ims','_coloc_spots_data.csv');
%     filename2=regexprep(filename2,'.czi','_coloc_spots_data.csv');
%     filename2=regexprep(filename2,'.lsm','_coloc_spots_data.csv');
%     filename2=regexprep(filename2,'.tif','_coloc_spots_data.csv');
%     filename2=regexprep(filename2,'.tiff','_coloc_spots_data.csv');
%     if size(filename)==size(filename2)
%         if filename2==filename
%             [filename2 pathname2]=uiputfile('*.csv','could not replace extension, please choose a location and name for .csv');
%         end    
%     end 
%     comma=',';
%     fid = fopen(filename2, 'w');
%     fprintf(fid, '%s', filename);
%     fprintf(fid, '%s', comma);
% 
%     fprintf(fid, '%g', Nspots(1));
%     fprintf(fid, '%s', comma);
% 
%     fprintf(fid, '%g', Nspots(2));
%     fprintf(fid, '%s', comma);
% 
%     fprintf(fid, '%g', Nspots(3));
%     fprintf(fid, '%s', comma);
% 
%     fprintf(fid, '%g', n1Colocs2);
%     fprintf(fid, '%s', comma);
% 
%     fprintf(fid, '%g', n2Colocs3);
%     fprintf(fid, '%s', comma);
% 
%     fprintf(fid, '%g', n1Colocs3);
%     fprintf(fid, '%s', comma);  
% 
%     fprintf(fid, '%g', n12Colocs23);
%     fprintf(fid, '%s', comma);
% 
%     fprintf(fid, '%s', [spotName{1},'-',spotName{2},'-',spotName{3}]);
%     fprintf(fid, '%s', comma);
% 
%     fprintf(fid, '%s\n', ' ');
% 
% 
%     fclose(fid);
%     for i=1:nSpotsObj
%         newFileName{i}=[spotName{i} '.csv'];
%     end    
%     filename3=regexprep(filename2,'.csv',newFileName{1});
%     filename4=regexprep(filename2,'.csv',newFileName{2});
%     filename5=regexprep(filename2,'.csv',newFileName{3});
%     if nSpotsObj>3
%         filename6=regexprep(filename2,'.csv',newFileName{4});
%     end
%     if nSpotsObj>4
%         filename7=regexprep(filename2,'.csv',newFileName{5});
%     end
%     if nSpotsObj>5
%         filename8=regexprep(filename2,'.csv',newFileName{6});
%     end    
%     
% end
end

