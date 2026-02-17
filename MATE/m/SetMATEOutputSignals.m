function SetMATEOutputSignals(sub,MATE,sensorV,sensorI,switchnames)

%sub='RouteTestModel_12_a/InternalSourceMixer1';

L=find_system(sub,'FindAll','on','type','line');

delete_line(L);
%nb_state nb_input nb_output nb_switches nb_nodes nb_permutations nb_ItypePorts

Block1 = 'Demux'; % Replace with your block name
Block2 = 'BusCreator'; % Replace with your block name
Block3 = 'SigOut'; % Replace with your block name
Block4 = 'Yout'; % Replace with your block name
%%MATE.sizes: nb_state nb_input nb_output nb_switches nb_nodes
nb_Vsen=size(sensorV,1);
nb_Isen=size(sensorI,1);
nbsensor=MATE.sizes(3)-MATE.sizes(4)-MATE.sizes(5);

ind=1;
%set_param([sub '/Demux'], 'Outputs',['[' num2str(MATE.sizes(4)) ' ' str ' ' num2str(MATE.sizes(5)) ']' ])
set_param([sub '/Demux'], 'Outputs', num2str(MATE.sizes(3)) );
set_param([sub '/BusCreator'], 'Inputs', num2str(MATE.sizes(3)) );

for i=1:MATE.sizes(4)
    h=add_line(sub, [Block1 '/' num2str(ind)], [Block2 '/' num2str(ind)]);
    w=get_param(switchnames{i},'Name');

    set_param(h,'Name',['SW_voltage:_' num2str(ind) '_' w]); 
    ind=ind+1;
end
for i=1:nb_Vsen
    h=add_line(sub, [Block1 '/' num2str(ind)], [Block2 '/' num2str(ind)]);
    try
        nn=get_param( get_param(sensorV{i,1},'Parent'),'Parent');
    catch
        nn=sensorV{i,1};
    end
    w=get_param(nn,'Name');

    set_param(h,'Name',['Voltage:_' num2str(ind) '_' w]); 
    ind=ind+1;
end
for i=1:nb_Isen
    h=add_line(sub, [Block1 '/' num2str(ind)], [Block2 '/' num2str(ind)]);
    try
        nn=get_param( get_param(sensorI{i,1},'Parent'),'Parent');
    catch
        nn=sensorI{i,1};
    end
    w=get_param(nn,'Name');

    set_param(h,'Name',['current:_' num2str(ind) '_' w]); 
    ind=ind+1;
end
for i=1:MATE.sizes(5)
    h=add_line(sub, [Block1 '/' num2str(ind)], [Block2 '/' num2str(ind)]);
    set_param(h,'Name',['MATEOut:_' num2str(ind) '_internal']); 
    ind=ind+1;
end
% reconnect output of Bus Creator
add_line(sub, [Block2 '/1'], [Block3 '/1']);
add_line(sub, [Block4 '/1'], [Block1 '/1']);




