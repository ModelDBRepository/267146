function [syn_tab_asc syn_tab_desc axon_tab_asc axon_tab_desc] = generate_axons_and_make_synapses_L(cell_types, rc, dorsal_dendrite, ventral_dendrite)
%
global cell_colours;
global side_shift;
global total_number_of_cells;
%
%%
syn_tab_asc=[];
syn_tab_desc=[];
axon_tab_asc=[];
axon_tab_desc=[];
    %
    asc=0;
    desc=1;
    als=2; % axon on left side
    ars=3; % axon on right side
    %
    tt(1:8000)=0;
%
 for i=1:total_number_of_cells
     cell_type = cell_types(i);
    %
        %do_ascending = true; % false; 
    % if not dIN or mn
     switch cell_type
% %         case 1
% %             [initial_depth,india] = get_initial_axon_pos_L(idhistogram_rb_asc, idcentre_rb_asc); % desc same as asc for these, doesn't matter
% % % ASC %%%%%           
% %             [coord_asc, synapse_indices_asc, synapse_depths_asc, synapse_xs_asc] = make_axon_L(initial_depth, india, rc(i), RB_params_asc, iahistogram_rb_asc, iacentre_rb_asc, idhistogram_rb_asc, idcentre_rb_asc, lhistog_rb_asc, lhistcentre_rb_asc, iaangles_rb_asc, dorsal_dendrite, ventral_dendrite, true, cell_types, i,rc); % true/is ascending                                    
% %             tt(1:8000)=0;
% %             tt1=[i, asc, cell_types(i), als, length(coord_asc), coord_asc];
% %             tt(1:length(tt1))=tt1;
% %             axon_tab_asc=[axon_tab_asc; tt];
% %             syn_tab_asc= [syn_tab_asc;[i*ones(1,length(synapse_indices_asc))' cell_types(i)*ones(1,length(synapse_indices_asc))' synapse_indices_asc' cell_types(synapse_indices_asc)' synapse_xs_asc' synapse_depths_asc']];  
% % % DESC %%%%%            
% %             [coord_desc, synapse_indices_desc, synapse_depths_desc, synapse_xs_desc] = make_axon_L(initial_depth, india, rc(i), RB_params_desc, iahistogram_rb_desc, iacentre_rb_desc, idhistogram_rb_desc, idcentre_rb_desc, lhistog_rb_desc, lhistcentre_rb_desc, iaangles_rb_desc, dorsal_dendrite, ventral_dendrite, false, cell_types, i,rc);
% %             tt(1:8000)=0;
% %             tt1=[i, desc, cell_types(i), als, length(coord_desc), coord_desc];
% %             tt(1:length(tt1))=tt1;
% %             axon_tab_desc=[axon_tab_desc; tt];
% %             syn_tab_desc= [syn_tab_desc;[i*ones(1,length(synapse_indices_desc))' cell_types(i)*ones(1,length(synapse_indices_desc))' synapse_indices_desc' cell_types(synapse_indices_desc)' synapse_xs_desc' synapse_depths_desc']];
            case 1
            [coord_asc, synapse_indices_asc, synapse_depths_asc, synapse_xs_asc, coord_desc, synapse_indices_desc, synapse_depths_desc, synapse_xs_desc] = Prototype_make_axon_rb_L(rc(i), dorsal_dendrite, ventral_dendrite, cell_types, i,rc); % true/is ascending   
            tt(1:8000)=0;
            tt1=[(i), asc, cell_types(i), als, length(coord_asc), coord_asc];
            tt(1:length(tt1))=tt1;
            axon_tab_asc=[axon_tab_asc; tt];
            if length(synapse_indices_asc)>0
              syn_tab_asc= [syn_tab_asc;[(i)*ones(1,length(synapse_indices_asc))' cell_types(i)*ones(1,length(synapse_indices_asc))' synapse_indices_asc' cell_types(synapse_indices_asc)' synapse_xs_asc' synapse_depths_asc']];     
            end;
            tt(1:8000)=0;
            tt1=[(i),desc, cell_types(i), als, length(coord_desc), coord_desc];
            tt(1:length(tt1))=tt1;
            axon_tab_desc=[axon_tab_desc; tt];
            if length(synapse_indices_desc)>0
              syn_tab_desc= [syn_tab_desc;[(i)*ones(1,length(synapse_indices_desc))' cell_types(i)*ones(1,length(synapse_indices_desc))' synapse_indices_desc' cell_types(synapse_indices_desc)' synapse_xs_desc' synapse_depths_desc']];
            end;
% % %%%%%%%%%%%%            
% %             %rectangle('Position',[rc(i)-0.1,initial_depth-0.1, 0.2, 0.2],'FaceColor',cell_colours(cell_types(i),1:3),'Curvature',[1,1])
% %         case 2
% %             [initial_depth,india] = get_initial_axon_pos_L(idhistogram_dlc_asc, idcentre_dlc_asc); % desc same as asc for these, doesn't matter
% % % ASC %%%%%%%            
% %             [coord_asc, synapse_indices_asc, synapse_depths_asc, synapse_xs_asc] = make_axon_L_op(initial_depth, india, rc(i), dlc_params_asc, iahistogram_dlc_asc, iacentre_dlc_asc, idhistogram_dlc_asc, idcentre_dlc_asc, lhistog_dlc_asc, lhistcentre_dlc_asc, iaangles_dlc_asc, dorsal_dendrite, ventral_dendrite, true, cell_types, i,rc); % true/is ascending            
% %             tt(1:8000)=0;
% %             tt1=[i, asc, cell_types(i), ars, length(coord_asc), coord_asc];
% %             tt(1:length(tt1))=tt1;
% %             axon_tab_asc=[axon_tab_asc; tt];
% %             syn_tab_asc= [syn_tab_asc;[(i)*ones(1,length(synapse_indices_asc))' cell_types(i)*ones(1,length(synapse_indices_asc))' (synapse_indices_asc+side_shift)' cell_types(synapse_indices_asc+side_shift)' synapse_xs_asc' synapse_depths_asc']];  
% % % Desc %%%%%%            
% %             [coord_desc, synapse_indices_desc, synapse_depths_desc, synapse_xs_desc] = make_axon_L_op(initial_depth, india, rc(i), dlc_params_desc, iahistogram_dlc_desc, iacentre_dlc_desc, idhistogram_dlc_desc, idcentre_dlc_desc, lhistog_dlc_desc, lhistcentre_dlc_desc, iaangles_dlc_desc, dorsal_dendrite, ventral_dendrite, false, cell_types, i,rc);
% %             tt(1:8000)=0;
% %             tt1=[i, desc, cell_types(i), ars, length(coord_desc), coord_desc];
% %             tt(1:length(tt1))=tt1;
% %             axon_tab_desc=[axon_tab_desc; tt];
% %             syn_tab_desc= [syn_tab_desc;[(i*ones(1,length(synapse_indices_desc)))' (cell_types(i)*ones(1,length(synapse_indices_desc)))'  (synapse_indices_desc+side_shift)' cell_types(synapse_indices_desc+side_shift)' synapse_xs_desc' synapse_depths_desc']];
% % %%%%%%%%%%%%            
            case 2
            [coord_asc, synapse_indices_asc, synapse_depths_asc, synapse_xs_asc, coord_desc, synapse_indices_desc, synapse_depths_desc, synapse_xs_desc] = Prototype_make_axon_dlc_L(rc(i), dorsal_dendrite, ventral_dendrite, cell_types, i,rc); % true/is ascending   
            tt(1:8000)=0;
            tt1=[(i), asc, cell_types(i), als, length(coord_asc), coord_asc];
            tt(1:length(tt1))=tt1;
            axon_tab_asc=[axon_tab_asc; tt];
            if length(synapse_indices_asc)>0
              syn_tab_asc= [syn_tab_asc;[(i)*ones(1,length(synapse_indices_asc))' cell_types(i)*ones(1,length(synapse_indices_asc))' (synapse_indices_asc+side_shift)' cell_types(synapse_indices_asc+side_shift)' synapse_xs_asc' synapse_depths_asc']];     
            end;
            tt(1:8000)=0;
            tt1=[(i),desc, cell_types(i), als, length(coord_desc), coord_desc];
            tt(1:length(tt1))=tt1;
            axon_tab_desc=[axon_tab_desc; tt];
            if length(synapse_indices_desc)>0
              syn_tab_desc= [syn_tab_desc;[(i)*ones(1,length(synapse_indices_desc))' cell_types(i)*ones(1,length(synapse_indices_desc))' (synapse_indices_desc+side_shift)' cell_types(synapse_indices_desc+side_shift)' synapse_xs_desc' synapse_depths_desc']];
            end;
% %
% %             %rectangle('Position',[rc(i)-0.1,initial_depth-0.1+100, 0.2, 0.2],'FaceColor',cell_colours(cell_types(i),1:3),'Curvature',[1,1])
% %         case 3
% %             [initial_depth,india] = get_initial_axon_pos_L(idhistogram_ain_asc, idcentre_ain_asc); % desc same as asc for these, doesn't matter
% % % ASC %%%%%%%            
% %             [coord_asc, synapse_indices_asc, synapse_depths_asc, synapse_xs_asc] = make_axon_L(initial_depth, india, rc(i), aIN_params_asc, iahistogram_ain_asc, iacentre_ain_asc, idhistogram_ain_asc, idcentre_ain_asc, lhistog_ain_asc, lhistcentre_ain_asc, iaangles_ain_asc, dorsal_dendrite, ventral_dendrite, true, cell_types, i,rc); % true/is ascending            
% %             tt(1:8000)=0;
% %             tt1=[i, asc, cell_types(i), als, length(coord_asc), coord_asc];
% %             tt(1:length(tt1))=tt1;
% %             axon_tab_asc=[axon_tab_asc; tt];
% %             syn_tab_asc= [syn_tab_asc;[i*ones(1,length(synapse_indices_asc))' cell_types(i)*ones(1,length(synapse_indices_asc))' synapse_indices_asc' cell_types(synapse_indices_asc)' synapse_xs_asc' synapse_depths_asc']];
% % % DSC %%%%%%            
% %             [coord_desc, synapse_indices_desc, synapse_depths_desc, synapse_xs_desc] = make_axon_L(initial_depth, india, rc(i), aIN_params_desc, iahistogram_ain_desc, iacentre_ain_desc, idhistogram_ain_desc, idcentre_ain_desc, lhistog_ain_desc, lhistcentre_ain_desc, iaangles_ain_desc, dorsal_dendrite, ventral_dendrite, false, cell_types, i,rc);
% %             tt(1:8000)=0;
% %             tt1=[i, desc, cell_types(i), als, length(coord_desc), coord_desc];
% %             tt(1:length(tt1))=tt1;
% %             axon_tab_desc=[axon_tab_desc; tt];
% %             syn_tab_desc= [syn_tab_desc;[i*ones(1,length(synapse_indices_desc))' cell_types(i)*ones(1,length(synapse_indices_desc))' synapse_indices_desc' cell_types(synapse_indices_desc)' synapse_xs_desc' synapse_depths_desc']];
            case 3
            [coord_asc, synapse_indices_asc, synapse_depths_asc, synapse_xs_asc, coord_desc, synapse_indices_desc, synapse_depths_desc, synapse_xs_desc] = Prototype_make_axon_aIN_L(rc(i), dorsal_dendrite, ventral_dendrite, cell_types, i,rc); % true/is ascending   
            tt(1:8000)=0;
            tt1=[(i), asc, cell_types(i), als, length(coord_asc), coord_asc];
            tt(1:length(tt1))=tt1;
            axon_tab_asc=[axon_tab_asc; tt];
            if length(synapse_indices_asc)>0
              syn_tab_asc= [syn_tab_asc;[(i)*ones(1,length(synapse_indices_asc))' cell_types(i)*ones(1,length(synapse_indices_asc))' synapse_indices_asc' cell_types(synapse_indices_asc)' synapse_xs_asc' synapse_depths_asc']];     
            end;
            tt(1:8000)=0;
            tt1=[(i),desc, cell_types(i), als, length(coord_desc), coord_desc];
            tt(1:length(tt1))=tt1;
            axon_tab_desc=[axon_tab_desc; tt];
            if length(synapse_indices_desc)>0
              syn_tab_desc= [syn_tab_desc;[(i)*ones(1,length(synapse_indices_desc))' cell_types(i)*ones(1,length(synapse_indices_desc))' synapse_indices_desc' cell_types(synapse_indices_desc)' synapse_xs_desc' synapse_depths_desc']];
            end;
% % %%%%%%%%%%%%            
% %             %rectangle('Position',[rc(i)-0.1,initial_depth-0.1, 0.2, 0.2],'FaceColor',cell_colours(cell_types(i),1:3),'Curvature',[1,1])
% %        case 4
% %             [initial_depth,india] = get_initial_axon_pos_L(idhistogram_cin_asc, idcentre_cin_asc); % desc same as asc for these, doesn't matter
% % % ASC %%%%%%%            
% %             [coord_asc, synapse_indices_asc, synapse_depths_asc, synapse_xs_asc] = make_axon_L_op(initial_depth, india, rc(i), cIN_params_asc, iahistogram_cin_asc, iacentre_cin_asc, idhistogram_cin_asc, idcentre_cin_asc, lhistog_cin_asc, lhistcentre_cin_asc, iaangles_cin_asc, dorsal_dendrite, ventral_dendrite, true, cell_types, i,rc); % true/is ascending            
% %             tt(1:8000)=0;
% %             tt1=[i, asc, cell_types(i), ars, length(coord_asc), coord_asc];
% %             tt(1:length(tt1))=tt1;
% %             axon_tab_asc=[axon_tab_asc; tt];
% %             syn_tab_asc= [syn_tab_asc;[(i)*ones(1,length(synapse_indices_asc))' cell_types(i)*ones(1,length(synapse_indices_asc))' (synapse_indices_asc+side_shift)' cell_types(synapse_indices_asc+side_shift)' synapse_xs_asc' synapse_depths_asc']];
% % % Desc %%%%%%            
% %             [coord_desc, synapse_indices_desc, synapse_depths_desc, synapse_xs_desc] = make_axon_L_op(initial_depth, india, rc(i), cIN_params_desc, iahistogram_cin_desc, iacentre_cin_desc, idhistogram_cin_desc, idcentre_cin_desc, lhistog_cin_desc, lhistcentre_cin_desc, iaangles_cin_desc, dorsal_dendrite, ventral_dendrite, false, cell_types, i,rc);
% %             tt(1:8000)=0;
% %             tt1=[i, desc, cell_types(i), ars, length(coord_desc), coord_desc];
% %             tt(1:length(tt1))=tt1;
% %             axon_tab_desc=[axon_tab_desc; tt];
% %             syn_tab_desc= [syn_tab_desc;[(i)*ones(1,length(synapse_indices_desc))' cell_types(i)*ones(1,length(synapse_indices_desc))' (synapse_indices_desc+side_shift)' cell_types(synapse_indices_desc+side_shift)' synapse_xs_desc' synapse_depths_desc']];
% % %%%%%%%%%%%%            
% %             %rectangle('Position',[rc(i)-0.1,initial_depth-0.1+100, 0.2, 0.2],'FaceColor',cell_colours(cell_types(i),1:3),'Curvature',[1,1])
            case 4
            [coord_asc, synapse_indices_asc, synapse_depths_asc, synapse_xs_asc, coord_desc, synapse_indices_desc, synapse_depths_desc, synapse_xs_desc] = Prototype_make_axon_cIN_L(rc(i), dorsal_dendrite, ventral_dendrite, cell_types, i,rc); % true/is ascending   
            tt(1:8000)=0;
            tt1=[(i), asc, cell_types(i), als, length(coord_asc), coord_asc];
            tt(1:length(tt1))=tt1;
            axon_tab_asc=[axon_tab_asc; tt];
            if length(synapse_indices_asc)>0
              syn_tab_asc= [syn_tab_asc;[(i)*ones(1,length(synapse_indices_asc))' cell_types(i)*ones(1,length(synapse_indices_asc))' (synapse_indices_asc+side_shift)' cell_types(synapse_indices_asc+side_shift)' synapse_xs_asc' synapse_depths_asc']];     
            end;
            tt(1:8000)=0;
            tt1=[(i),desc, cell_types(i), als, length(coord_desc), coord_desc];
            tt(1:length(tt1))=tt1;
            axon_tab_desc=[axon_tab_desc; tt];
            if length(synapse_indices_desc)>0
              syn_tab_desc= [syn_tab_desc;[(i)*ones(1,length(synapse_indices_desc))' cell_types(i)*ones(1,length(synapse_indices_desc))' (synapse_indices_desc+side_shift)' cell_types(synapse_indices_desc+side_shift)' synapse_xs_desc' synapse_depths_desc']];
            end;
%
           case 5
           flag_asc=0;
                if rc(i) <=850
                  work_space=rand;
                  if work_space <0.85
                     flag_asc=1;
                  end;
                else
                   if (rc(i) >850 && rc(i)<=1400)
                      work_space=rand;
                      if work_space <0.6 
                         flag_asc=1;
                      end;
                   end;
                end;
            [coord_asc, synapse_indices_asc, synapse_depths_asc, synapse_xs_asc, coord_desc, synapse_indices_desc, synapse_depths_desc, synapse_xs_desc] = Prototype_make_axon_dIN_L(rc(i), dorsal_dendrite, ventral_dendrite, cell_types, i,rc,flag_asc); % true/is ascending   
            %[coord_desc, synapse_indices_desc, synapse_depths_desc, synapse_xs_desc] = make_axon_L(initial_depth, india, rc(i), dIN_params, iahistogram_din, iacentre_din, idhistogram_din, idcentre_din, lhistog_din, lhistcentre_din, iaangles_din, dorsal_dendrite, ventral_dendrite, false, cell_types, i,rc);   
            tt(1:8000)=0;
            tt1=[i,desc, cell_types(i), als, length(coord_desc), coord_desc];
            tt(1:length(tt1))=tt1;
            axon_tab_desc=[axon_tab_desc; tt];
            if length(synapse_indices_desc)>0
              syn_tab_desc= [syn_tab_desc;[i*ones(1,length(synapse_indices_desc))' cell_types(i)*ones(1,length(synapse_indices_desc))' synapse_indices_desc' cell_types(synapse_indices_desc)' synapse_xs_desc' synapse_depths_desc']];
            end;  
                if flag_asc==1
                  tt(1:8000)=0;
                  tt1=[i, asc, cell_types(i), als, length(coord_asc), coord_asc];
                  tt(1:length(tt1))=tt1;
                  axon_tab_asc=[axon_tab_asc; tt];
                  if length(synapse_indices_asc)>0
                     syn_tab_asc= [syn_tab_asc;[i*ones(1,length(synapse_indices_asc))' cell_types(i)*ones(1,length(synapse_indices_asc))' synapse_indices_asc' cell_types(synapse_indices_asc)' synapse_xs_asc' synapse_depths_asc']];     
                  end;
                end; 
 %   end;
            %
% %             case 6
% %             [initial_depth,india] = get_initial_axon_pos_L(idhistogram_mn, idcentre_mn); % desc same as asc for these, doesn't matter            
% %             %synapse_indices_asc = [];
% %             %synapse_xs_asc=[];
% %             %synapse_depths_asc=[];
% % % THERE ARE No ASC SO Desc only
% %             [coord_desc, synapse_indices_desc, synapse_depths_desc, synapse_xs_desc] = make_axon_L(initial_depth, india, rc(i), mn_params, iahistogram_mn, iacentre_mn, idhistogram_mn, idcentre_mn, lhistog_mn, lhistcentre_mn, iaangles_mn, dorsal_dendrite, ventral_dendrite, false, cell_types, i,rc);
% %             tt(1:8000)=0;
% %             tt1=[i,desc, cell_types(i), als, length(coord_desc), coord_desc];
% %             tt(1:length(tt1))=tt1;
% %             axon_tab_desc=[axon_tab_desc; tt];
% %             %syn_tab_asc= [syn_tab_asc;[i*ones(1,length(synapse_indices_asc))' cell_types(i)*ones(1,length(synapse_indices_asc))' synapse_indices_asc' cell_types(synapse_indices_asc)' synapse_xs_asc' synapse_depths_asc']];
% %             syn_tab_desc= [syn_tab_desc;[i*ones(1,length(synapse_indices_desc))' cell_types(i)*ones(1,length(synapse_indices_desc))' synapse_indices_desc' cell_types(synapse_indices_desc)' synapse_xs_desc' synapse_depths_desc']];    
% % %%%%%%%%%%%%            
% %             %rectangle('Position',[rc(i)-0.1,initial_depth-0.1, 0.2, 0.2],'FaceColor',cell_colours(cell_types(i),1:3),'Curvature',[1,1])
            case 6
            [coord_asc, synapse_indices_asc, synapse_depths_asc, synapse_xs_asc, coord_desc, synapse_indices_desc, synapse_depths_desc, synapse_xs_desc] = Prototype_make_axon_mn_L(rc(i), dorsal_dendrite, ventral_dendrite, cell_types, i,rc); % true/is ascending   
            tt(1:8000)=0;
            tt1=[(i), asc, cell_types(i), als, length(coord_asc), coord_asc];
            tt(1:length(tt1))=tt1;
            axon_tab_asc=[axon_tab_asc; tt];
            if length(synapse_indices_asc)>0
            end;
            tt(1:8000)=0;
            tt1=[(i),desc, cell_types(i), als, length(coord_desc), coord_desc];
            tt(1:length(tt1))=tt1;
            axon_tab_desc=[axon_tab_desc; tt];
            if length(synapse_indices_desc)>0
              syn_tab_desc= [syn_tab_desc;[(i)*ones(1,length(synapse_indices_desc))' cell_types(i)*ones(1,length(synapse_indices_desc))' synapse_indices_desc' cell_types(synapse_indices_desc)' synapse_xs_desc' synapse_depths_desc']];
            end;
%    end
        case 7
            [coord_asc, synapse_indices_asc, synapse_depths_asc, synapse_xs_asc, coord_desc, synapse_indices_desc, synapse_depths_desc, synapse_xs_desc] = Prototype_make_axon_dla_L(rc(i), dorsal_dendrite, ventral_dendrite, cell_types, i,rc);
            %[initial_depth,india] = get_initial_axon_pos_L(idhistogram_dla, idcentre_dla); % desc same as asc for these, doesn't matter            
%             synapse_indices_desc = [];
%             synapse_xs_desc=[];
%             synapse_depths_desc=[];
% THERE ARE No Desc SO ASC only
            %[coord_asc, synapse_indices_asc, synapse_depths_asc, synapse_xs_asc] = make_axon_L(initial_depth, india, rc(i), dla_params, iahistogram_dla, iacentre_dla, idhistogram_dla, idcentre_dla, lhistog_dla, lhistcentre_dla, iaangles_dla, dorsal_dendrite, ventral_dendrite, true, cell_types, i,rc); % true/is ascending            
            tt(1:8000)=0;
            tt1=[i, asc, cell_types(i), als, length(coord_asc), coord_asc];
            tt(1:length(tt1))=tt1;
            axon_tab_asc=[axon_tab_asc; tt];
            if length(synapse_indices_asc)>0
              syn_tab_asc= [syn_tab_asc;[i*ones(1,length(synapse_indices_asc))' cell_types(i)*ones(1,length(synapse_indices_asc))' synapse_indices_asc' cell_types(synapse_indices_asc)' synapse_xs_asc' synapse_depths_asc']];
            end;
            tt(1:8000)=0;
            tt1=[(i),desc, cell_types(i), als, length(coord_desc), coord_desc];
            tt(1:length(tt1))=tt1;
            axon_tab_desc=[axon_tab_desc; tt];
            if length(synapse_indices_desc)>0
              syn_tab_desc= [syn_tab_desc;[(i)*ones(1,length(synapse_indices_desc))' cell_types(i)*ones(1,length(synapse_indices_desc))' synapse_indices_desc' cell_types(synapse_indices_desc)' synapse_xs_desc' synapse_depths_desc']];
            end;
            %%%%%%%%%%%%            
            %rectangle('Position',[rc(i)-0.1,initial_depth-0.1, 0.2, 0.2],'FaceColor',cell_colours(cell_types(i),1:3),'Curvature',[1,1])
         %otherwise
     end
     %
end     

%end
%% Truncate axon tables (delete zeros)
 % if cell_types ~= 6 
% 
   axon_asc_L_nmax=max(axon_tab_asc(:,5))+5;
   axon_tab_asc=axon_tab_asc(:,1:axon_asc_L_nmax);
%  end
%
   axon_desc_L_nmax=max(axon_tab_desc(:,5))+5;
   axon_tab_desc=axon_tab_desc(:,1:axon_desc_L_nmax);
%
%
end
