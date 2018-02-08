function [LORALL, PTf] = iniPT(f_min, f_max, y, c)

% 
maxlevel = 10;

PTf = []; %zeros(cut_num, 8);

% first cut point
newitem(1) = f_min;
newitem(2) = f_max;
newitem(3) = (f_max-f_min)/2 + f_min; 

newitem(4) = countnum(y, newitem(1), newitem(3)); % n0_1
newitem(5) = countnum(y, newitem(3), newitem(2)); % n0_2
newitem(6) = newitem(4); % n1_1
newitem(7) = newitem(5); % n1_2

newitem(8) = evaluteLOR(newitem(4:7), c, 1); % LOR
LORALL     = newitem(8);

newitem(9:10) = 0;
newitem(11) = c;

PTf        = [PTf; newitem];
queuesize  = 1;

% width-first tree traversal 
intervel  = [newitem(1) newitem(3); newitem(3) newitem(2)];
parentidx = [1 1];

newitem_leftchild  = [];
newitem_rightchild = [];

for lev = 2 : maxlevel
    
    for cut_num_level = 1 : 2^(lev-1-1)
        
        % left child
        newitem_leftchild(1) = intervel(1, 1); % f_min
        newitem_leftchild(2) = intervel(1, 2); % f_max
        newitem_leftchild(3) = (newitem_leftchild(2)-newitem_leftchild(1))/2 + newitem_leftchild(1); % f_cut
        
        newitem_leftchild(4) = countnum(y, newitem_leftchild(1), newitem_leftchild(3)); % n0_1
        newitem_leftchild(5) = countnum(y, newitem_leftchild(3), newitem_leftchild(2)); % n0_2
        newitem_leftchild(6) = newitem_leftchild(4); % n1_1
        newitem_leftchild(7) = newitem_leftchild(5); % n1_2
        
        newitem_leftchild(8) = evaluteLOR(newitem_leftchild(4:7), c, lev); % LOR
        LORALL               = LORALL + newitem_leftchild(8);
        
        newitem_leftchild(9:10) = 0;        
        newitem_leftchild(11)   = c * lev^2;
        
        % right child
        newitem_rightchild(1) = intervel(2, 1);
        newitem_rightchild(2) = intervel(2, 2);
        newitem_rightchild(3) = (newitem_rightchild(2)-newitem_rightchild(1))/2 + newitem_rightchild(1);
        
        newitem_rightchild(4) = countnum(y, newitem_rightchild(1), newitem_rightchild(3)); % n0_1
        newitem_rightchild(5) = countnum(y, newitem_rightchild(3), newitem_rightchild(2)); % n0_2
        newitem_rightchild(6) = newitem_rightchild(4); % n1_1
        newitem_rightchild(7) = newitem_rightchild(5); % n1_2
        
        newitem_rightchild(8) = evaluteLOR(newitem_rightchild(4:7), c, lev);
        LORALL                = LORALL + newitem_rightchild(8);
        
        newitem_rightchild(9:10) = 0;
        newitem_rightchild(11)   = c * lev^2;
        
        % add to queue
        PTf(parentidx(1), 9)  = queuesize + 1; % left child
        PTf(parentidx(2), 10) = queuesize + 2; % right child
        
        parentidx(1:2)   = []; % remove the first two
        parentidx        = [parentidx queuesize + 1 queuesize + 1 queuesize + 2 queuesize + 2];
        
        intervel(1:2, :) = []; % remove the first two
        intervel         = [intervel; [newitem_leftchild(1) newitem_leftchild(3)]; [newitem_leftchild(3) newitem_leftchild(2)]];
        intervel         = [intervel; [newitem_rightchild(1) newitem_rightchild(3)]; [newitem_rightchild(3) newitem_rightchild(2)]];
        
        PTf       = [PTf; newitem_leftchild; newitem_rightchild];
        queuesize = queuesize + 2;
    end
    
end

% output  

end


function num = countnum(list, f_min, f_max)
    num = numel(find(list >= f_min & list < f_max));    
end

function LOR = evaluteLOR(num, c, level)
    alpha0 = c * level^2;
    alpha1 = alpha0;
    
    contrib_num = - betaln(alpha0,alpha1) + betaln(alpha0 + sum(num(1:2)), alpha1 + sum(num(3:4)));
    contrib_den = - betaln(alpha0,alpha1) - betaln(alpha0,alpha1)...
                    + betaln(alpha0 + num(1), alpha1 + num(3)) + betaln(alpha0 + num(2), alpha1 + num(4));
    
    LOR = contrib_num - contrib_den;
end


