function [MATE, MATERT,blk]= MakeMate2(modele, disc)

%Ts=25e-6;
SILENT=0;  % 0 = less disp()
opt=0;

% 23 juillet  works on RTScape_netlist23juillet.slx
%zip('backup.zip', {'*.m', '*.mat'});


model_top=get_param(modele,'Parent');

if isempty(model_top)
    Ts=evalin('base', get_param(bdroot,'FixedStep'));
else
    Ts=evalin('base', get_param(model_top,'FixedStep'));
end

%% DISCRETISE
if 1==0
disc=1; %1 backeuler  2 trap  5 pade(2,3)   3: pade(2,3_+BE_input
end

%modele=bdroot;
disp(['Running node analysis on ' modele]);

% at the start comment both MATE source
mateblk=[find_system(modele, 'regexp', 'on', 'FollowLinks', 'on','LookUnderMasks','all','MaskType','MATE Node')];
matetmp={};
for i=1:size(mateblk,1)
    try 
        a=get_param(mateblk{i},'conty');
        % if does not have this property, this is a single pphase block
    catch
        x=size(matetmp,1);
        matetmp{x+1,1}= mateblk{i};
    end
end
mateblk_sld=matetmp;   % a patch until we modify all test case with SLD MATE blocks with conty parameter

mateblk_1ph=[find_system(modele, 'regexp', 'on', 'FollowLinks', 'on','LookUnderMasks','all','MaskType','MATE Node','conty', 'Single-Phase')];



MakeMATEInterface([],mateblk_sld,mateblk_1ph,[],[],[],[],[],[],'commentMATE');

blk=find_system(modele, 'regexp', 'on', 'FollowLinks', 'on','LookUnderMasks','all','ReferenceBlock','ee_lib');
blk=[blk;find_system(modele, 'regexp', 'on', 'FollowLinks', 'on','LookUnderMasks','all','ReferenceBlock','fl_lib')]; % base ssc lib block

nb_blk=size(blk,1);

% easy way to find MATE type from the first block encountered...
%[mateblk] = MakeMATEInterface(blk,mateblk,[],[],[],'step1');





han=zeros(nb_blk,3);
rlc=zeros(1,6);
src=zeros(1,6);
sensors=zeros(1,6);



blkmask={
'Electrical Reference',[1],[1,0],-1,{'LConn1'} ; 
 ['Grounded Neutral' 10 '(Three-Phase)'],[1 0],[1,0],-1,{'TBD'};
'Resistor',[1 2],[1 1],-1,{'LConn1', 'RConn1'};                              % block type, electric ports indice for LConn and RConn in this order
'Inductor',[1 2],[1 1],-1,{'TBD'};                              % followed by [number of LConn, number of RConn]
'Capacitor',[1 2],[1 1],-1,{'TBD'};                             % followed by modeling option if any (-1 if none)
%'RLC (Three-Phase)',[1 2 3 4 5 6],[3 3],'ee.passive.rlc_assemblies.rlc.Xabc'; % ABC ports  % modelingOption = get_param(gcbh, 'ComponentPath') 'ee.passive.rlc_assemblies.rlc.abc' get_param(gcb, 'port_option'):'ee.enum.threePhasePort.expanded'
'RLC (Three-Phase)',[1 2 3 4 5 6],[3 3],'ee.enum.threePhasePort.expanded',{'TBD'}; % ABC ports  % modelingOption = get_param(gcbh, 'ComponentPath') 'ee.passive.rlc_assemblies.rlc.abc' get_param(gcb, 'port_option'):'ee.enum.threePhasePort.expanded'
'RLC (Three-Phase)',[1 2],[1 1],        'ee.passive.rlc_assemblies.rlc.abc',{'TBD'};  % sld ports 
'AC Current Source',[1 2],[1 1],-1,{'TBD'};
'AC Voltage Source',[1 2],[1 1],-1,{'TBD'};
'DC Voltage Source',[1 2],[1 1],-1,{'TBD'};   %10
 ['Voltage' 10 'Source' 10 '(Three-Phase)'],[1 2],[1 1],'ee.sources.voltage.abc',{'TBD'};  %SLD
 'Voltage Sensor',[1 3],[1 1],-1,{'LConn1', 'RConn2'};
 'Current Sensor',[1 3],[1 1],-1,{'LConn1', 'RConn2'}; 
 [    'Current and Voltage' 10 'Sensor (Three-Phase)'],[1 4],[1 3],'ee.sensors.vi_sensor.abc' ,{'TBD'}; % SLD %this block has 2 Rconn for outputs meas
 'Phase Splitter',[1 2 3 4],[1,3],-1,{'TBD'};  % 4 conns
 ['Circuit Breaker' 10 '(Three-Phase)'],[2 3 4 5 6 7],[3 3],'ee.switches.circuit_breaker.ps.Xabc',{'TBD'};
 ['Current' 10 'Source' 10 '(Three-Phase)'],[1 2],[1 1],'ee.sources.current.abc',{'TBD'};  %SLD
 ['Coupled Lines' 10 '(Three-Phase)'],[1 2 3 4 5 6],[3 3], 'ee.passive.lines.coupled_lines.Xabc',{'TBD'};  % mutual inductance
 ['Nonlinear' 10 'Transformer'],[1 2 3 4],[2 2], -1,{'TBD'};   %#19 linear xfo   %21 aout : switch position with mutual inductance
 'Switch',[1 3],[1 1],-1,{'TBD'};  %20
 'Voltage Source',[1 2],[1 1],-1,{'LConn1','RConn1'};
 'Current Source',[1 2],[1 1],-1,{'TBD'};
  ['Circuit Breaker' 10 '(Three-Phase)'],[2 3],[1 1],'ee.switches.circuit_breaker.ps.abc',{'TBD'}; %23
  ['Two-Winding' 10 'Transformer' 10 '(Three-Phase)'],[1 2],    [1 1], 'ee.passive.transformers.two_winding_transformer.abc',{'TBD'};  %24  NOTE: Y port not available yet
  ['Controlled Current' 10 'Source'],[1 3],[1 1],-1,{'LConn1', 'RConn2'};  %25
};


SLD_block_index=...
    {[7 11 14 17 23 24 ], ...   %{1} SLD blocks  
    [6 16 18],...        %{2} 3-ph expended blocks  28 aout add 18 mutual inductance!
    [15],...          %{3} phase splitter
    [13 14],...        %{4} current sensors
    [3 4 5 6 7 18 19],...     %{5} rlc branches
    [11],...            % 6: source with RL SLD
    [12],...              %7 voltage sensors
    [-1],...       % 8group of same type, different port option,
    [8 22 25]};    % current source only

nb_blktype=size(blkmask,1);





% find the type of block
qq=[];
for i=1:nb_blk
    flag=0;
    for j=1:nb_blktype
        if isequal(get_param(blk{i,1},'MaskType'),blkmask{j,1})
            if isequal(blkmask{j,4},-1)
                qq=j;
                flag=1;
                break;
            else                
                % cond=isequal(get_param(blk{i,1}, 'ComponentPath'),blkmask{j,4});
                % if cond==0
                % try
                %    cond=isequal(get_param(blk{i,1}, 'port_option'),blkmask{j,4});
                % catch
                % end
                % end
                    
                cond1=0; cond2=0;
                try 
                    cond1=isequal(get_param(blk{i,1}, 'ComponentPath'),blkmask{j,4}) ;                          
                catch
                   cond1=0;
                end
                try
                    cond2=isequal(get_param(blk{i,1}, 'port_option'),blkmask{j,4}) ; 
                catch
                    cond2=0;
                end
                    
                %if isequal(get_param(blk{i,1}, 'ComponentPath'),blkmask{j,4})
                if cond1==1 | cond2==1
                
                    if isempty(find(SLD_block_index{8}==j))
                            qq=j;
                            flag=1;
                            break;
                    else
                        tmp=get_param(blk{i,1},'Ports');
                        if isequal(tmp(6:7),blkmask{j,3})
                            qq=j;
                            flag=1;
                            break;
                        end
                    end
                end
            end
        end
    end
    if flag==0
        error(['blocktype is not supported yet:' blk{i,1}]);
    end
    blk{i,2}=qq;
    for k=1: length(find(blkmask{qq,2}>0))  % noeud PMIO >0
        blk{i,2+k}=-1;   % initial node number attribution to -1
    end
    if isempty(blk{i,2})
        error('blocktype is null');
    end

end

% reorder to put grounds first
blk1={}; i=1;
blk2={}; j=1;
blk3={}; jj=1;
blk4={}; jjj=1;
cols=size(blk,2);
for k=1:nb_blk
    if (blk{k,2}==1) || (blk{k,2}==2) % ground &  % 3-phase ground
        for m=1:cols
            blk1{i,m}=blk{k,m};
        end
        i=i+1;
    elseif blk{k,2}==19
        for m=1:cols
            blk3{jj,m}=blk{k,m};
        end
        jj=jj+1;
    elseif blk{k,2}==18  %mutual inductance
        for m=1:cols
            blk4{jjj,m}=blk{k,m};
        end
        jjj=jjj+1;

    else
        for m=1:cols
            blk2{j,m}=blk{k,m};
        end
        j=j+1;
    end
end
blk=[blk1 ;blk2; blk3; blk4];
blkgnd=blk1;
%blk=[blk1 ; blk3; blk2];


net=0;
nets=[];
all_block={};
all_types=cell(nb_blk,2);
physicalBlockList={};
PortTypeList={};

for i=1:nb_blk
    %if isequal(blk{i,1},'MATE_HDLsolver_inverter_1/MATEREF/Voltage Source')
    if isequal(blk{i,1}, 'MATE_HDLsolver_inverter_1/MATEREF/Voltage Source')
        sss=1
    end
    if i>=6
        sss=1
    end

    for j=1:size(blkmask{blk{i,2},5},2)
        contyp=blkmask{blk{i,2},5}{j};
        disp([blk{i,1} '  ' contyp]);




        [physicalBlockList_tmp, PortTypeList_tmp]=findConnectedBlocks(blk{i,1},contyp);
        physicalBlockList_tmp=[physicalBlockList_tmp;blk{i,1}];
        PortTypeList_tmp=[ PortTypeList_tmp ; contyp];
        gnd_blk=0;


        flag=0;
        for k=1:length(nets)
            for m=1:size(physicalBlockList,1)
                for n=1:size(physicalBlockList{k},1)
                    for o=1:size(physicalBlockList_tmp,1)
                        if isequal(physicalBlockList_tmp{o},physicalBlockList{k}{n})
                            if isequal(PortTypeList_tmp{o},PortTypeList{k}{n})
                                flag=1;
                                break;
                            end
                        end
                    end
                    if flag==1
                        break;
                    end
                end
                if flag==1
                    break;
                end
            end
            if flag==1
                break;
            end
        end
        if flag==0
            net=net+1;
            physicalBlockList{net}=physicalBlockList_tmp;
            PortTypeList{net}=PortTypeList_tmp;
            if i<=size(blkgnd,1)
                nets=[nets; 0];
            else
                nets=[nets; net];
            end
        end
    end
end

        

%% assign the node in the blk list
for i=1:length(nets)
    if i==7
        aaa=1;
    end
    for j=1:size(physicalBlockList{i},1)
        flag=0;
        for k=1:nb_blk
            thisBlock=blk{k,1};
            if isequal(thisBlock,'MATE_HDLsolver_inverter_1/MATEREF/resload1')
                xxx=1;
            end
            thisBlockType=blk{k,2};
            thisBlock_NbConn=size(blkmask{thisBlockType,5},2);
            if isequal(thisBlock,physicalBlockList{i}{j})
                for m=1:thisBlock_NbConn
                    if isequal(blkmask{thisBlockType,5}{m},PortTypeList{i}{j})
                        blk{k,2+m}=nets(i);
                        flag=1;
                        break;
                    end
                end
            end
            if flag==1
                break;
            end
        end
    end
end
MATE={}; MATERT={};



NodeNumber=max(nets)+1;





%if i>nb_blk
    disp('BEFORE Expension ');
    % add rl of 3 ph source here.
    %TOBE REFINED 
    blk = AddRLsrcNodes(blk,blkmask,SLD_block_index, NodeNumber+4, blkgnd); 
    break_main=1;
    netlistprint(blk);
    disp(' ');
    disp('AFTER SLD Expension ')
    newblk = ExpendSLDNodes(blk,blkmask,SLD_block_index, 1000);  %expand 3-phase blocks
    blk=newblk;
    netlistprint(blk);
    disp(' ');
    disp('AFTER Splitter Expension ')
    blk=MergeSpliterNode(blk,blkmask,SLD_block_index, 1000);
    netlistprint(blk)
    disp(' ');
    disp('AFTER Current sensor node merger ')
    %[blk, Vnode,cur_br_all, sensorV,sensorI] = MergeCurrentSensorNodes(blk,SLD_block_index);
    [blk, Vnode,cur_br_all, sensorV,sensorI,Vmeas,VSensorNames] = MergeCurrentSensorNodes(blk,SLD_block_index);
    disp('Current sensors and blocks ')
    if SILENT
    for j=1:size(cur_br_all,1)
        disp(['sensor: ' cur_br_all{j,1}]);   
        for k=1:length(cur_br_all{j,2})
            disp(['       right: ' num2str(j) '  ' blk{cur_br_all{j,2}(k),1}]);
        end
        for k=1:length(cur_br_all{j,3})
            disp(['       left: ' num2str(j) '  ' blk{cur_br_all{j,3}(k),1}]);
        end
    end
    end
        netlistprint(blk)
        disp(' ');
        disp('Make SPS netlist ')
        [rlc, source, switches, sensors, rlcnames,sourcenames,source_parameter,switchnames,switchtype,switchVf,switch_parameter] = MakeSPSNetlist(blk,blkmask,1000);
        %netlistprint(blk)
        %[yout,y_type,output_sig] = BuildSPSYout(Vnode, switches, cur_br_all,rlcnames,sourcenames,switchnames);
        % if opt(1)==1
        %     %add here the L anc C node voltage
        %     Vnode=[];
        %     for k=1:size(rlc,1)
        %         if rlc(k,5)>0
        %             if rlc(k,4)>0 || rlc(k,6)>0
        %                 error('Under FPGA option (opt(1)=1)  R L C branches can be mixed (ex: RL)')
        %             end
        %             Vnode=[Vnode; rlc(k,1) rlc(k,2)];
        %         end
        %     end
        % end


        [yout,y_type,output_sig,allSensorNames] = BuildSPSYout(Vnode, switches, cur_br_all,rlcnames,sourcenames,switchnames,Vmeas,VSensorNames);
        %y_type=[zeros(1,size(switches,1)) zeros(1,length(Vnode))];  % %%%  added Vnode manually!
        disp('Build MATE interface ')
        [mateblk,source,yout,y_type,sourcenames,nb_nodal_nodes,portType,src_order] = MakeMATEInterface(blk,mateblk_sld,mateblk_1ph,rlc,source,switches,yout,y_type,sourcenames,'findMATEtype');

        disp('Compute ABCD matrices ')
        unit='OMU'; % Units = ohms, mH, and uF
        % %add part of switches on top of source
        % if ~isempty(switches)
        %     swsrc=[switches(:,1:2) ones(size(switches,1),1) zeros(size(switches,1),3)];
        %     source=[swsrc; source];
        % end

        % Open file that contains power_statespace output information
        %fid=fopen('power_circ2ss.net','w');
        %rlc=[rlc(77:end,:); rlc(1:76,:)] ;
         % tmp=1:size(rlc,1);
         % rlc=[rlc tmp'];
         % tmp=ones(size(source,1),1)+1;
         % source=[source tmp];
 
        % rlc=[rlc(end-2:end,:);rlc(1:end-3,:)];
        %[A,B,C,D,states,x0,x0sw,rlsw,u,x,y,freq,Asw,Bsw,Csw,Dsw,Hlin]=power_statespace(rlc,switches,source,[],yout,y_type,unit);
        %[A,B,C,D,states,x0,x0sw,rlsw,u,x,y,freq,Asw,Bsw,Csw,Dsw,Hlin]=mate_statespace2(rlc,switches,source,[],yout,y_type,unit);
        [A,B,C,D,states,x0,x0sw,rlsw,u,x,y,freq,Asw,Bsw,Csw,Dsw,Hlin,Outputs]=mate_statespace3(rlc,switches,source,[],yout,y_type,unit);
        %Ts=50e-6; should be in the parameter input
        disp('Discretize things... ')
        %nb_nodal_nodes=3;
        %portType=source(end-nb_nodal_nodes+1:end,3);

        %% DISCRETISE METHOD
        if ~isempty(switches)
            [Adp, B1dp, Cdp, Ddp, Yp, Cinj_p, Dinj_p, z_p, x_p, B2dp]= MATESwitchPermute(A,B,C,D,switches(:,4)',Ts,nb_nodal_nodes,portType,disc);
        else
            [Adp, B1dp, Cdp, Ddp, Yp, Cinj_p, Dinj_p, z_p, x_p, B2dp]= MATESwitchPermute(A,B,C,D,[],Ts,nb_nodal_nodes,portType,disc);
        end
        statenames=[];   %for later
        u0=zeros(size(source,1),1);
        x0=zeros(size(Adp{1},1),1);  % for later, if required
        disp('Prepare S-function data ')
        MATE=MATE_Sfun_DataFormat(rlc, source, switches, yout,Adp,B1dp,B2dp,Cdp,Ddp, Yp, Cinj_p, Dinj_p, z_p, x_p,mateblk,u0,x0, statenames, rlcnames, sourcenames,switchtype,nb_nodal_nodes,portType, Ts,switchVf);
        MATERT=MATE_Mex_DataFormat(A,B,C,D, switches,switchtype,disc,u0,x0, Ts,MATE.sizes,switchVf);
        disp('Reorder the MATE source from external connexions ')
        try
        RouteMATESources(modele, src_order,'MATE_Equivalent/InputElectricMixer');  %% must change the name if in a subsystem
        catch
            warning('No rewrite order of source done!!!')
        end
        disp('Adding Simulink source in place of Simscape ones ')
        try
        [srcParam,src_Selector]=AddSources(source_parameter,size(switches,1),nb_nodal_nodes);  %srcParam is use directly as a parameter in the model
        catch
            warning('No adds of source done!!!')
        end
        MATE.srcParam=srcParam;
        MATE.src_Selector=src_Selector;
        MATE.source_parameter=source_parameter;
        MATE.switch_parameter=switch_parameter;
        MATE.src_order=src_order;
        MATE.model=modele;
        MATE.Outputs=Outputs;
        MATE.allSensorNames=allSensorNames;

        % if opt(1)==1
        %     disp('Generating MATE FPGA fixed-admittance matrices')
        %     MATE=MateConvertForFPGA(MATE);
        %     disp('FPGA routine: finished')
        %     return;
        % end

        disp('Routing Simulink switch signal in place of Simscape ones ') 
        try
        [swParam, masterSwSelOrder]=AddSwitches(modele,switchnames,'InternalSwitchMixer');
        catch
            warning('Error in switch routing... continuing!!!')
        end
        disp('Disabling MATE sources in order to enable simulation from external network')
        MakeMATEInterface([],mateblk_sld,mateblk_1ph,[],[],[],[],[],[],'commentMATE');
        try
            disp('Routing MATE Outputs')   
            SetMATEOutputSignals([model_top '/MATE_Equivalent/MATEOutputSelector'],MATE,sensorV,sensorI,switchnames);
        catch
            warning('Error in output routing... continuing!!!')
        end
    %end


    % print out

end



