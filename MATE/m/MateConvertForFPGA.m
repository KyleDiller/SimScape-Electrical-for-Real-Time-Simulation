function[MATE_SOC]=MateConvertForFPGA(MATE)
% format parameter for FPGA simulation

MAX_INPUTS=16;
MAX_OUTPUTS=16;
MAX_SWITCHES=8;

sourcenames=MATE.sourcenames(:,1);
MATE_SOC.MATE=MATE;

%find switch pairs

Elements=find_system(MATE.model, 'FollowLinks', 'on','LookUnderMasks','all','MaskType','Pejovic Element for FPGA simulation');
Inverters2L=find_system(MATE.model, 'FollowLinks', 'on','LookUnderMasks','all','MaskType','MATE Compensated Pejovic 2-level inverter');
Switches=find_system(MATE.model, 'FollowLinks', 'on','LookUnderMasks','all','MaskType','Pejovic Element for FPGA simulation','typ','Switch');
Inductors=find_system(MATE.model, 'FollowLinks', 'on','LookUnderMasks','all','MaskType','Pejovic Element for FPGA simulation','typ','Inductor');
Capacitors=find_system(MATE.model, 'FollowLinks', 'on','LookUnderMasks','all','MaskType','Pejovic Element for FPGA simulation','typ','Capacitor');
SWF_2Linverters=find_system(MATE.model, 'FollowLinks', 'on','LookUnderMasks','all','MaskType','MATE High impedance capable Switching Function inverter (2-level IGBT/Diode)');

order=[];
ordertype=[];
for i=1:size(SWF_2Linverters,1)
    vsen=[SWF_2Linverters{i} '/Vsensor'];
    isen=[SWF_2Linverters{i}  '/Isensor'];
    vsrc=[SWF_2Linverters{i}  '/CVS'];
    isrc=[SWF_2Linverters{i}  '/CCS'];
    flag=0;
    for j=1:size(MATE.sourcenames,1)
        if isequal(vsrc,MATE.sourcenames{j,1})
            flag=1;
            break;
        end

    end
    if flag==0
        error('MateConvertForFPGA:30 v')
    end
    order=[order j];
    ordertype=[ordertype 2];
    flag=0;
    for j=1:size(MATE.sourcenames,1)
        if isequal(isrc,MATE.sourcenames{j,1})
            flag=1;
            break;
        end

    end
    if flag==0
        error('MateConvertForFPGA:30 i')
    end
    order=[order j];
    ordertype=[ordertype 2 3];
end

for i=1:size(Inverters2L,1)
    upp=[Inverters2L{i} '/upper/CCS1'];
    loo=[Inverters2L{i} '/lower/CCS1'];
    orr=[];
    orrs=[0 0];
    for j=1:size(MATE.sourcenames,1)
        if isequal(upp,MATE.sourcenames{j,1})
            orr(1)=j;
            orrs(1)=1;
        end
        if isequal(loo,MATE.sourcenames{j,1})
            orr(2)=j;
            orrs(2)=1;
        end
        if sum(orrs)==2
            break;
        end
    end
    if sum(orrs)<2
        error('MateConvertForFPGA:1')
    end
    order=[order orr];
    ordertype=[ordertype 1 1];
end
for i=1:size(Switches,1)
    upp=[Switches{i} '/CCS1'];
    flag=0;
    for j=1:size(MATE.sourcenames,1)
        if isequal(upp,MATE.sourcenames{j,1})
            flag=1;
            break;
        end

    end
    if flag==0
        error('MateConvertForFPGA:4')
    end
    order=[order j];
    ordertype=[ordertype 1];
end
for i=1:size(Inductors,1)
    upp=[Inductors{i} '/CCS1'];
    flag=0;
    for j=1:size(MATE.sourcenames,1)
        if isequal(upp,MATE.sourcenames{j,1})
            flag=1;
            break;
        end

    end
    if flag==0
        error('MateConvertForFPGA:2')
    end
    order=[order j];
    ordertype=[ordertype 1];

end
for i=1:size(Capacitors,1)
    upp=[Capacitors{i} '/CCS1'];
    flag=0;
    for j=1:size(MATE.sourcenames,1)
        if isequal(upp,MATE.sourcenames{j,1})
            flag=1;
            break;
        end

    end
    if flag==0
        error('MateConvertForFPGA:3')
    end
    order=[order j];
    ordertype=[ordertype 1];

end

% find the rest of source and put it at the end of the inputs.
notPejoInputs=[];
for i=1:size(MATE.source,1)
    if ~any(order==i, "all")
        notPejoInputs=[notPejoInputs i];
    end
end



outorder=[];

for i=1:size(SWF_2Linverters,1)
    vsen=[SWF_2Linverters{i}  '/Vsensor'];
    isen=[SWF_2Linverters{i}  '/Isensor'];
    flag=0;
    for j=1:size(MATE.allSensorNames,1)
        if isequal(vsen,MATE.allSensorNames{j,1})
            flag=1;
            break;
        end

    end
    if flag==0
        error('MateConvertForFPGA:31 vsensor')
    end
    outorder=[outorder j];
    flag=0;
    for j=1:size(MATE.allSensorNames,1)
        if isequal(isen,MATE.allSensorNames{j,1})
            flag=1;
            break;
        end
    end
    if flag==0
        error('MateConvertForFPGA:31 isensor')
    end
    outorder=[outorder j];

end
for i=1:length(order)
    weHaveMatch=0;
    if ordertype(i)==1
        nodes=MATE.source(order(i),1:2);
        % find the index in yout
        for j=1:size(MATE.yout,1)
            isU_n=strfind(MATE.yout(j,:),'U_n');
            if ~isempty(isU_n)
                loc=strfind(MATE.yout(j,:),'_');
                outnodes=[str2num(MATE.yout(j,loc(2)+1:end)) str2num(MATE.yout(j,loc(1)+2:loc(2)-1))];  % attention: order is reversed in yout from source
            end
            if length(outnodes)==2
                if nodes(1)==outnodes(1)  && nodes(2)==outnodes(2)
                    weHaveMatch=1;
                    outorder=[outorder j];
                    break;
                end
            end
        end
        if weHaveMatch==0
            error('MateConvertForFPGA:5')
        end
    end
end

% find the rest of yout and put it at the end of the outputs
notPejoOutputs=[];
for i=1:size(MATE.yout,1)
    if ~any(outorder==i, "all")
        notPejoOutputs=[notPejoOutputs i];
    end
end

% permute Inputs of D

order=[order notPejoInputs];
outorder=[outorder notPejoOutputs];
DD=MATE.Ddp{1}(outorder,order);
MATE_SOC.new_sourcenames=sourcenames(order);
MATE_SOC.new_outputnames=MATE.allSensorNames(outorder);
[a,b]=size(DD);
DD=[DD zeros(a,MAX_INPUTS-b)];
DD=[DD; zeros(MAX_INPUTS-a, MAX_INPUTS)];
MATE_SOC.D=DD;
MATE_SOC.orderin=order;
MATE_SOC.orderout=outorder;







