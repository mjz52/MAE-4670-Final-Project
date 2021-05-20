% Michael Zakoworotny
% J integral function
% Uses a direct integration method along the sides of a rectangular contour
% in a rectangular mesh

% cont_size - contour size, number of elements per half-side of countour
% path

% d_ind - load step increment index to compute the J integral for

% Currently supports cracks on the lower face of a rectangular domain 
% (face 1), and only Q4 elements

function J = J_integral(meshStruct, boundStruct, globalSystem, cont_size, d_ind)

    nCoords = meshStruct.nCoords;
    elCon = meshStruct.elCon;
    nnpe = meshStruct.nnpe;
    gatherMat = meshStruct.gatherMat;
    glU = globalSystem.u_k(:,d_ind);

    % Get stress evaluated at each node
    [~, stress_field] = calNodalStrainStress(glU,meshStruct);
    
    % Get list of elements and face numbers in integral
    % Node num of crack tip
    tip_num = boundStruct.nodes(end).Nodes(1); % this will be true on all faces
    [el_surr,~] = find(meshStruct.elCon == tip_num); % elements containing crack node
    el_surr = sort(el_surr);
    numElx = length(boundStruct.elements(3).Elems);
    numEly = length(boundStruct.elements(2).Elems);
    
    contour_list = [];
    face_list = [];
    % Contour 1 - right
    el_cont_1 = (el_surr(2)-1+cont_size) : numElx : (el_surr(2)-1+cont_size)+numElx*(cont_size-1);
    contour_list = [contour_list, el_cont_1];
    face_list = [face_list, repelem(2,1,length(el_cont_1))];
    
    % Contour 2 - upper
    el_cont_2 = flip((el_surr(1)+1-cont_size)+numElx*(cont_size-1) : el_cont_1(end));
    contour_list = [contour_list, el_cont_2];
    face_list = [face_list, repelem(3,1,length(el_cont_2))];
    
    % Contour 3 - left
    el_cont_3 = flip(el_surr(1)+1-cont_size : numElx : el_cont_2(end));
    contour_list = [contour_list, el_cont_3];
    face_list = [face_list, repelem(4,1,length(el_cont_3))];
    
    % Compute J integral for each contour side
    J = 0;
    for e = 1:length(contour_list)
        elNum = contour_list(e);
        faceNum = face_list(e);
        % Get unit normal, node id's edge, and nodal coords in parent
        % domain
        switch faceNum
            case 2
                n = [1 0]';
                node1 = elCon(elNum,2);
                node2 = elCon(elNum,3);
                qp1 = [1,-1];
                qp2 = [1,1];
            case 3
                n = [0 1]';
                node1 = elCon(elNum,3);
                node2 = elCon(elNum,4);
                qp1 = [1,1];
                qp2 = [-1,1];
            case 4
                n = [-1 0]';
                node1 = elCon(elNum,4);
                node2 = elCon(elNum,1);
                qp1 = [-1,1];
                qp2 = [-1,-1];
        end
        
        % Elemental information
        lcU = glU(gatherMat(elNum,:));
        u_elem = [lcU(1:2:end),lcU(2:2:end)];
        xynode = nCoords(elCon(elNum,:),:); % nodal coordinates - [x, y]
        
        % Length of element side
        el_side = sqrt((nCoords(node2,1)-nCoords(node1,1))^2 + (nCoords(node2,2)-nCoords(node1,2))^2);
        
        % First point
        dNdXiEta = dNmatrix(qp1,nnpe);
        JofXiEta_t = dNdXiEta*xynode; % Jacobian of element evaluated at point
        dNdXY = JofXiEta_t\dNdXiEta;
        dNdx1 = dNdXY(1,:);
        % Get stretch tensor
        xy_def = xynode + u_elem;
        Fe = xy_def'*dNdXY';
        Ce = Fe'*Fe ;
        Ee = 1/2*(Ce - eye(2));
        % Get strain energy density
        [~,~,W] = meshStruct.Material.constLaw(Ee, meshStruct.Material, e);

        sig = stress_field(node1,:);
        sig_tens = [sig(1), sig(3); sig(3), sig(2)];
        t = sig_tens*n; % traction vector
        J_i = W*n(1) - t' * sum(dNdx1.*u_elem', 2);
        
        % Second point
        dNdXiEta = dNmatrix(qp2,nnpe);
        JofXiEta_t = dNdXiEta*xynode; % Jacobian of element evaluated at point
        dNdXY = JofXiEta_t\dNdXiEta;
        dNdx1 = dNdXY(1,:);
        % Get stretch tensor
        xy_def = xynode + u_elem;
        Fe = xy_def'*dNdXY';
        Ce = Fe'*Fe ;
        Ee = 1/2*(Ce - eye(2));
        % Get strain energy density
        [~,~,W] = meshStruct.Material.constLaw(Ee, meshStruct.Material, e);

        sig = stress_field(node2,:);
        sig_tens = [sig(1), sig(3); sig(3), sig(2)];
        t = sig_tens*n; % traction vector
        % Evaluate integrand at node
        J_ip1 = W*n(1) - t' * sum(dNdx1.*u_elem', 2);
        % Add contribution to sum using area integration
        J = J + 1/2*(J_i + J_ip1)*el_side;
    end
    
    J = J*2; % multiply by two because of symmetry

end