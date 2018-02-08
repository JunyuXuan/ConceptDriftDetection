function [h, post, stats, PTf] = incPTtest(PTf, olddatum, newdatum, c, threshold)

% tree width-first traversal
nodequeue = 1;

% checkstatus = 2 -- both data need to be checked; 
% checkstatus = 1 -- only check olddatum; 
% checkstatus = 0 -- only check newdatum; 
checkstatus = 2; 

while ~isempty(nodequeue)
    
    current_node_idx = nodequeue(1); 
    
    if current_node_idx == 0
        nodequeue(1)   = [];
        checkstatus(1) = [];
        continue;
    end
    
    checkchild       = checkstatus(1);
    node_info        = PTf(current_node_idx, :);
    cut_value        = node_info(3);
    alpha            = node_info(11);
    
    switch checkchild
        case 0
            
            %disp('check newdatum') 
            if newdatum >= cut_value
                updateLOR   = node_info(8) + log(alpha + node_info(5) ) - log(2*alpha + sum(node_info(4:5))) ...
                                  + log(2*alpha + sum(node_info(6:7)) ) - log(alpha + node_info(7)); 
                              
                PTf(current_node_idx, 8) = updateLOR;
                PTf(current_node_idx, 7) = node_info(7) + 1;
                
                nodequeue   = [nodequeue node_info(10)];
                checkstatus = [checkstatus 0];                
            else
                updateLOR   = node_info(8) + log(alpha + node_info(4) ) - log(2*alpha + sum(node_info(4:5))) ...
                                  + log(2*alpha + sum(node_info(6:7)) ) - log(alpha + node_info(6)); 
                              
                PTf(current_node_idx, 8) = updateLOR;
                PTf(current_node_idx, 6) = node_info(6) + 1;
                
                nodequeue   = [nodequeue node_info(9)];
                checkstatus = [checkstatus 0];                
            end
            
        case 1
            
            %disp('check olddatum')
            if  olddatum >= cut_value
                updateLOR   = node_info(8) + log(2*alpha + sum(node_info(4:5)) - 1) - log(alpha + node_info(5) - 1) ...
                                  + log(alpha + node_info(7) - 1) - log(2*alpha + sum(node_info(6:7)) - 1) ; 
                
                PTf(current_node_idx, 8) = updateLOR;
                PTf(current_node_idx, 7) = node_info(7) - 1;
                                              
                nodequeue   = [nodequeue node_info(10)];
                checkstatus = [checkstatus 1];                
            else
                updateLOR   = node_info(8) + log(2*alpha + sum(node_info(4:5)) - 1) - log(alpha + node_info(4) - 1) ...
                                  + log(alpha + node_info(6) - 1) - log(2*alpha + sum(node_info(6:7)) - 1) ; 
                
                PTf(current_node_idx, 8) = updateLOR;
                PTf(current_node_idx, 6) = node_info(6) - 1;
                
                nodequeue   = [nodequeue node_info(9)];
                checkstatus = [checkstatus 1];
            end
            
        case 2
            
            %disp('check two nodes')
            if olddatum < cut_value && newdatum < cut_value
                
                nodequeue   = [nodequeue node_info(9)];
                checkstatus = [checkstatus 2];
            else
                if olddatum >= cut_value && newdatum >= cut_value                    
                    nodequeue   = [nodequeue node_info(10)];
                    checkstatus = [checkstatus 2];
                else
                    nodequeue = [nodequeue node_info(9) node_info(10)];
                    
                    if olddatum >= cut_value                        
                        updateLOR   = node_info(8) + log(alpha + node_info(4)) - log(alpha + node_info(6)) ...
                                            + log(alpha + node_info(7) - 1) - log(alpha + node_info(5) - 1) ; 
                        
                        PTf(current_node_idx, 8) = updateLOR;
                        PTf(current_node_idx, 7) = node_info(7) - 1;
                        PTf(current_node_idx, 6) = node_info(6) + 1;   
                                        
                        checkstatus = [checkstatus 0 1];
                    else
                        updateLOR   = node_info(8) + log(alpha + node_info(4)) - log(alpha + node_info(6)) ...
                                            + log(alpha + node_info(7) - 1) - log(alpha + node_info(5) - 1) ; 
                        
                        PTf(current_node_idx, 8) = updateLOR;
                        PTf(current_node_idx, 7) = node_info(7) + 1;
                        PTf(current_node_idx, 6) = node_info(6) - 1;  
                        
                        checkstatus = [checkstatus 1 0];
                    end
                end
            end
    end
    
    nodequeue(1)   = [];
    checkstatus(1) = [];    
end

% output    
LOR       = sum(PTf(:, 8));
post      = 1/(1+exp(-LOR));
h         = post < threshold;
stats.LOR = LOR;
stats.c   = c;


