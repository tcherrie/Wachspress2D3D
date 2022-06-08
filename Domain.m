
classdef Domain
    %Domain Class for domain objects
    %
    % To create a domain, use the constructor obj = Domain(type)
    %
    % For 0D, 1D and 2D regular polytopes, type is a number which
    % indicates the number of vertices :
    %  type = 1 : 1 vertex = a dot (0D)
    %  type = 2 : 2 vertices = a line (1D)
    %  type = n > 2 : a regular 2D polygon with n vertices
    %
    % For 3D polytopes, type is a string which denotes a specific
    % polytope. Available :
    %  "tetra" or "tetraedron" -> tetraedron
    %  "cube" -> cube
    %  "diamondN" -> diamond with a base of N vertices (N is an integer)
    %  "prismN" -> prism with a base of N vertices (N is an integer)
    %
    %Copyright (C) 2022 Théodore CHERRIERE (theodore.cherriere@ens-paris-saclay.fr)
    %This program is free software: you can redistribute it and/or modify
    %it under the terms of the GNU General Public License as published by
    %the Free Software Foundation, either version 3 of the License, or
    %any later version.
    %
    %This program is distributed in the hope that it will be useful,
    %but WITHOUT ANY WARRANTY; without even the implied warranty of
    %MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %GNU General Public License for more details.
    %
    %You should have received a copy of the GNU General Public License
    %along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
    properties
        vertices  % cartesian coordinates of vertices
        edges     % connectivity of edges in relation to vertices (3D domains only)
        facets    % connectivity of facets i.r.t edges(3D domains only)
        normals   % outgoing normal vectors i.r.t facets (3D domains only)
        vertices2facets   % connectivity of vertices i.r.t facets (3D domains only)
        edges2facets % connectivity of edges i.r.t vertices (3D domains only)
        dimension % dimension of the domain (0D, 1D, 2D or 3D)
        normal_fan % normal fan of the domain (computed automatically)
    end
    
    methods
        function obj = Domain(type)
            % obj = Domain(type)
            % Constructor of a domain
            % For 0D, 1D and 2D regular polytopes, type is a number which
            % indicates the number of vertices :
            %  type = 1 -> a dot (0D)
            %  type = 2 -> a line (1D)
            %  type = n > 2 -> a regular 2D polygon with n vertices
            %
            % For 3D polytopes, type is a string which denotes a specific
            % polytope. Available :
            %  "tetra" or "tetraedron" -> tetraedron
            %  "cube" -> cube
            %  "diamondN" -> diamond with a base of N vertices (N is an integer)
            %  "prismN" -> prism with a base of N vertices (N is an integer)
            % Feel free to add other types of 3D domains.
            
            if class(type)=="double" && type==1 % dimension 0 : a dot
                obj.vertices=0;
                obj.dimension=0;
                v=obj.vertices;
                
            elseif class(type)=="double" && type==2 % dimension 1 : a line
                obj.vertices=[-0.5;0.5];
                obj.dimension=1;
                v=obj.vertices;
                
            elseif class(type)=="double" % dimension 2 : a regular polygon
                if  length(type)==1
                    theta=linspace(0,2*pi,type+1)+(2*pi)/(2*type);
                    theta=shiftdim(theta(1:end-1),-2);
                    R=mult([cos(theta),-sin(theta);sin(theta),cos(theta)],[1;0]);
                    obj.vertices=reshape(R,[2,type]).';
                    v=obj.vertices;
                else
                    v=type;
                    obj.vertices=type;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% 3D
                
            elseif lower(type)=="tetra" || lower(type)=="tetraedron"
                theta=linspace(0,2*pi,3+1)+(2*pi)/(2*3);
                theta=shiftdim(theta(1:end-1),-2);
                R=mult([cos(theta),-sin(theta);sin(theta),cos(theta)],[1;0]);
                v1=reshape(R,[2,3]).';
                obj.vertices=[v1,zeros(3,1);[0 0 1.5]];
                obj.vertices=obj.vertices-mean(obj.vertices,1);
                v=obj.vertices;
                
                un(1,:)=cross(v(4,:)-v(1,:),v(2,:)-v(1,:)); un(2,:)=cross(v(4,:)-v(2,:),v(3,:)-v(2,:));
                un(3,:)=cross(v(4,:)-v(3,:),v(1,:)-v(3,:)); un(4,:)=cross(v(2,:)-v(1,:),v(3,:)-v(1,:));
                obj.normals=-un./vecnorm(un,2,2);
                
                obj.vertices2facets={[1,3,4];[2,1,4];[3,2,4];[1,2,3]};
                obj.edges=[1,2;2,3;3,1;4,1;4,2;4,3];
                obj.facets={[1,-5,4];[2,-6,5];[3,-4,6];[-1,-3,-2]};
                
                obj.edges2facets=zeros(length(obj.edges),2);
                for i=1:length(obj.facets)
                    edges=abs(obj.facets{i});
                    for j=1:length(edges)
                        if obj.edges2facets(edges(j),1)==0
                            obj.edges2facets(edges(j),1)=i;
                        else
                            obj.edges2facets(edges(j),2)=i;
                        end
                    end
                end
                
            elseif lower(type)=="cube"
                v=[0,0,0;1 0 0;0 1 0 ;1 1 0; 0 0 1; 1 0 1; 0 1 1; 1 1 1]; obj.vertices=v-mean(v,1);
                obj.vertices2facets={[1,4,6];[2,1,6];[4,3,6];[3,2,6];[5,4,1];[5,1,2];[5,3,4];[5,2,3]};
                
                un(1,:)=cross(v(2,:),v(5,:)); un(2,:)=cross(v(4,:)-v(2,:),v(6,:)-v(2,:));
                un(3,:)=cross(v(3,:)-v(4,:),v(8,:)-v(4,:)); un(4,:)=cross(v(1,:)-v(3,:),v(7,:)-v(3,:));
                un(5,:)=cross(v(6,:)-v(5,:),v(7,:)-v(5,:)); un(6,:)=cross(v(1,:)-v(2,:),v(4,:)-v(2,:));
                obj.normals=un./vecnorm(un,2,2);
                
                obj.edges=[1,2;2,4;4,3;3,1;5,6;6,8;8,7;7,5;1,5;2,6;4,8;3,7];
                obj.facets={[1,10,-5,-9];[2,11,-6,-10];[3,12,-7,-11];[4,9,-8,-12];...
                    [5,6,7,8];[1,2,3,4]};
                obj.edges2facets=zeros(length(obj.edges),2);
                for i=1:length(obj.facets)
                    edges=abs(obj.facets{i});
                    for j=1:length(edges)
                        if obj.edges2facets(edges(j),1)==0
                            obj.edges2facets(edges(j),1)=i;
                        else
                            obj.edges2facets(edges(j),2)=i;
                        end
                    end
                end
                v=obj.vertices;
                
            elseif contains(lower(type),"diamond")
                n  = str2double(extract(lower(type),regexpPattern("\d+")));
                theta=linspace(0,2*pi,n+1)+pi/n;
                theta=shiftdim(theta(1:end-1),-2);
                R=mult([cos(theta),-sin(theta);sin(theta),cos(theta)],[1;0]);
                v=reshape(R,[2,n]).'; v=[[v,zeros(n,1)];0,0,1;0,0,-1]; obj.vertices=v;
                for i=1:n
                    i1=i;
                    i2 = mod(i+n-2,n)+1;
                    i3 = i2+n;
                    i4 = i+n;
                    obj.vertices2facets{i} = [i1,i2,i3,i4];
                end
                obj.vertices2facets{n+1} = 1:n;
                obj.vertices2facets{n+2} = (2*n):-1:(n+1);
                for i=1:n
                    un(i,:)=cross(v(mod(i,n)+1,:)-v(i,:),v(n+1,:)-v(i,:));
                    un(i+n,:)=cross(v(n+2,:)-v(i,:),v(mod(i,n)+1,:)-v(i,:));
                end
                obj.normals=un./vecnorm(un,2,2);
                for i=1:n
                    obj.edges(i,:)=[i,mod(i,n)+1];
                    obj.edges(i+n,:)=[i,1+n];
                    obj.edges(i+2*n,:)=[2+n,i];
                end
                for i=1:n
                    i1 = i;
                    i2 = n + 1+ mod(i,n);
                    i3 = - n - i;
                    
                    i2b = -i2 - n;
                    i3b = -i3+n;
                    
                    obj.facets{i} = [i1,i2,i3];
                    obj.facets{i+n} = [i1,i2b,i3b];
                end
                obj.edges2facets=zeros(length(obj.edges),2);
                for i=1:length(obj.facets)
                    edges=abs(obj.facets{i});
                    for j=1:length(edges)
                        if obj.edges2facets(edges(j),1)==0
                            obj.edges2facets(edges(j),1)=i;
                        else
                            obj.edges2facets(edges(j),2)=i;
                        end
                    end
                end
                v=obj.vertices;
                
            elseif contains(lower(type),"prism")
                n  = str2double(extract(lower(type),regexpPattern("\d+")));
                theta=linspace(0,2*pi,n+1)+pi/n;
                theta=shiftdim(theta(1:end-1),-2);
                R=mult([cos(theta),-sin(theta);sin(theta),cos(theta)],[1;0]);
                v=reshape(R,[2,n]).'; v=[[v,-0.5*ones(n,1)];[v,0.5*ones(n,1)]]; obj.vertices=v;
                for i=1:n
                    i1=i;
                    i2 = mod(i+n-2,n)+1;
                    obj.vertices2facets{i} = [i1,i2,n+1];
                    obj.vertices2facets{i+n} = [i1,n+2,i2];
                end
                for i=1:n
                    i2 = mod(i,n)+1;
                    obj.edges(i,:) = [i,i2];
                    obj.edges(i+n,:) = [i+n,i2+n];
                    obj.edges(i+2*n,:)= [i,i+n];
                end
                for i=1:n
                    i1=i;
                    i2 = mod(i,n)+1+2*n;
                    i3 = -mod(i-1,n)-1-n;
                    i4 = -mod(i-1,n)-1-2*n;
                    obj.facets{i} = [i1,i2,i3,i4];
                end
                obj.facets{n+1}=1:n;
                obj.facets{n+2}=n+(1:n);
                for i=1:n
                    un(i,:)=cross(v(mod(i,n)+1,:)-v(i,:),v(n+i,:)-v(i,:));
                end
                un(n+1,:)=[0,0,-1];
                un(n+2,:)=[0,0,1];
                obj.normals=un./vecnorm(un,2,2);
                obj.edges2facets=zeros(length(obj.edges),2);
                for i=1:length(obj.facets)
                    edges=abs(obj.facets{i});
                    for j=1:length(edges)
                        if obj.edges2facets(edges(j),1)==0
                            obj.edges2facets(edges(j),1)=i;
                        else
                            obj.edges2facets(edges(j),2)=i;
                        end
                    end
                end
            end
            
            if size(obj.vertices,2)==2
                [cones_edge,cones_vertices,edges,vertices]=normalFan2D(v);
                obj.normal_fan.cones_edge=cones_edge;
                obj.normal_fan.cones_vertices=cones_vertices;
                obj.normal_fan.edges=edges;
                obj.normal_fan.vertices=vertices;
                obj.dimension=2;
                
                obj.edges=edges;
                normales=[0,1;-1 0]*(obj.vertices(obj.edges(:,2),:)-obj.vertices(obj.edges(:,1),:)).';
                obj.normals=((normales).')./vecnorm(normales.',2,2);
                
            elseif size(obj.vertices,2)==3
                [cones_facets,cones_edges,cones_vertices,facets,edges,vertices]=normalFan3D(obj);
                obj.normal_fan.cones_facets=cones_facets;
                obj.normal_fan.cones_edges=cones_edges;
                obj.normal_fan.cones_vertices=cones_vertices;
                obj.normal_fan.facets=facets;
                obj.normal_fan.edges=edges;
                obj.normal_fan.vertices=vertices;
                obj.dimension=3;
            end
        end
        
        %% Projection
        
        function val = projection(obj,val)
            % projection on the domain
            epsilon=1e-7;
            if obj.dimension==1
                val(val<-0.5)=-0.5;
                val(val>0.5)=0.5;
            elseif obj.dimension==2
                con_edg = obj.normal_fan.cones_edge;
                con_ver = obj.normal_fan.cones_vertices;
                edg = obj.normal_fan.edges;
                vert = obj.normal_fan.vertices;
                % Projection with the normal fan
                val=projection2D(val,con_ver,con_edg,vert,edg,epsilon);
                
            elseif obj.dimension==3
                con_fac = obj.normal_fan.cones_facets;
                con_edg = obj.normal_fan.cones_edges;
                con_ver = obj.normal_fan.cones_vertices;
                fac = obj.normal_fan.facets;
                edg = obj.normal_fan.edges;
                vert = obj.normal_fan.vertices;
                % Projection with the normal fan
                val=projection3D(val,obj.normals,con_fac,con_edg,con_ver,fac,edg,vert,epsilon);
            end
        end
        
        
        %% Display
        
        function plot_Wachspress(obj,n,legend_vertices)
            if nargin<=1 || isempty(n)
                n=1;
            end
            
            if nargin<=2 || isempty(legend_vertices)
                for i=1:length(obj.vertices)
                    legend_vertices(i)=convertCharsToStrings(strcat("v",num2str(i)));
                end
            end
            
            v_centre=obj.vertices-mean(obj.vertices,1);
            if size(v_centre,2)==1 && size(v_centre,1)==2 % 1D = 2 materials
                x=linspace(v_centre(1),v_centre(2),100).';
                w=wachspress(x,obj);
                w=w(n,:,:);
                plot(x,w(:));
                xlabel('x');
                ylabel('\omega(x)');
                for i=1:length(legend_vertices)
                    try
                        text(v_centre(i)*1.15,0,legend_vertices(i))
                    end
                end
            elseif size(v_centre,2)==2 %2D
                
                pgon=polyshape(v_centre(:,1)*0.99999,v_centre(:,2)*0.99999);
                [x,y]=meshgrid(linspace(-max(abs(v_centre(:))),max(abs(v_centre(:))),200),linspace(-max(abs(v_centre(:))),max(abs(v_centre(:))),200));
                in=pgon.isinterior(x(:),y(:));
                tr = delaunayTriangulation(x(in),y(in));
                w=wachspress([x(in),y(in)],obj);
                c=permute(w(n,:,:),[3 2 1]);
                trisurf(tr.ConnectivityList,tr.Points(:,1),tr.Points(:,2),c)
                colorbar off; maxval=max(tr.Points(:));
                axis([-maxval maxval -maxval maxval]*1.5);axis equal
                hold on; plot(pgon,'facealpha',0); grid on;
                shading interp
                for i=1:length(legend_vertices)
                    try
                        text(v_centre(i,1)*1.15,v_centre(i,2)*1.2,legend_vertices(i))
                    end
                end
                axis off
                hold off;
                light('Position',[-50 -15 10])
                material([0.8 0.5 0.3])
                lighting flat
                view(45,30)
                
            elseif   size(v_centre,2)==3 %3D
                
                shp=alphaShape(v_centre(:,1),v_centre(:,2),v_centre(:,3),Inf);
                [BF, P] = boundaryFacets(shp);
                [x,y]=meshgrid(0:0.01:1,0:0.01:1);
                tri1=delaunay(x,y);
                tri=[];
                nodes=[];
                for i=1:size(BF,1)
                    plan1=[0 0 0;1 0 0;0 1 0]+2;
                    plan2=[P(BF(i,:),1),P(BF(i,:),2),P(BF(i,:),3)];
                    M=get_transform_plan(plan1,plan2);
                    nodes_plan=M*[x(:).'+2;y(:).'+2;ones(1,length(x(:)))*2;ones(1,length(x(:)))];
                    tri=[tri;tri1+length(nodes)];
                    nodes=[nodes;nodes_plan(1:3,:).'];
                end
                nodes=nodes*0.99999;
                
                X0=sum([nodes(tri(:,1),1),nodes(tri(:,2),1),nodes(tri(:,3),1)],2)/3;
                Y0=sum([nodes(tri(:,1),2),nodes(tri(:,2),2),nodes(tri(:,3),2)],2)/3;
                Z0=sum([nodes(tri(:,1),3),nodes(tri(:,2),3),nodes(tri(:,3),3)],2)/3;
                in=shp.inShape(X0,Y0,Z0);
                tri=tri(in,:);
                w=wachspress(nodes,obj); w=w(n,:,:);
                w(w>1)=1; w(w<0)=0; % for points outside the polyhedron
                trisurf(tri,nodes(:,1),nodes(:,2),nodes(:,3),w(:))
                shading interp
                axis equal
                colorbar off; grid on;
                light('Position',[-50 -15 10])
                material([0.8 0.5 0.3])
                lighting flat
                view(45,30)
                axis off
                
                for i=1:length(legend_vertices)
                    try
                        text(v_centre(i,1)*1.15,v_centre(i,2)*1.2,v_centre(i,3)*1.2,legend_vertices(i))
                    end
                end
                
            end
            
            colormap(jet)
            colorbar
            
            if mod(n,10)==1
                legend(sprintf('%d^{st} basis function',n))
            elseif mod(n,10)==2
                legend(sprintf('%d^{nd} basis function',n))
            elseif mod(n,10)==3
                legend(sprintf('%d^{rd} basis function',n))
            else
                legend(sprintf('%d^{th} basis function',n))
            end
        end
        
        function plot(obj,legend_mat)
            if nargin<=1 || isempty(legend_mat)
                for i=1:length(obj.vertices)
                    legend_mat(i)=convertCharsToStrings(strcat("v",num2str(i)));
                end
            end
            
            couleur = legend2color(legend_mat);
            n_couleurs=length(obj.vertices);
            v_center=obj.vertices-mean(obj.vertices,1);
            
            if size(v_center,2)==1 && size(v_center,1)==2 % 1D = 2 materials
                x=linspace(v_center(1),v_center(2),100).';
                y=zeros(100,1) ;
                w=wachspress(x,obj);
                c=permute(sum(couleur(1:2,:).*(w.^5),1)+0.5,[3 2 1]);
                surf([x(:) x(:)], [y(:) y(:)], [x(:) x(:)],'FaceColor', 'none','EdgeColor', 'interp');
                colormap(c);
                axis(1.4*[-max(abs(v_center(:))),max(abs(v_center(:))),-max(abs(v_center(:))),max(abs(v_center(:)))])
                for i=1:length(legend_mat)
                    try
                        text(v_center(i)*1.15,0,legend_mat(i))
                    end
                end
                axis off
                
            elseif size(v_center,2)==2 %2D
                
                pgon=polyshape(v_center(:,1)*0.99999,v_center(:,2)*0.99999);
                [x,y]=meshgrid(linspace(-max(abs(v_center(:))),max(abs(v_center(:))),100),linspace(-max(abs(v_center(:))),max(abs(v_center(:))),100));
                in=pgon.isinterior(x(:),y(:));
                tr = delaunayTriangulation(x(in),y(in));
                nodes=tr.Points; elements=tr.ConnectivityList;
                X0=sum([nodes(elements(:,1),1),nodes(elements(:,2),1),nodes(elements(:,3),1)],2)/3;
                Y0=sum([nodes(elements(:,1),2),nodes(elements(:,2),2),nodes(elements(:,3),2)],2)/3;
                w=wachspress([X0,Y0],obj);
                
                c=permute(sum(couleur(1:n_couleurs,:).*(w.^5),1)+0.5,[3 2 1]);
                plotcolor2D(nodes,elements,(1:length(elements)));
                colormap(c); colorbar off; maxval=max(abs(nodes(:)));
                axis([-maxval maxval -maxval maxval]*1.5);axis equal
                hold on; plot(pgon,'facealpha',0); grid on;
                
                for i=1:length(legend_mat)
                    try
                        text(v_center(i,1)*1.15,v_center(i,2)*1.2,legend_mat(i))
                    end
                end
                axis off
                hold off;
                
            elseif size(v_center,2)==3 %3D
                
                shp=alphaShape(v_center(:,1),v_center(:,2),v_center(:,3),Inf);
                [BF, P] = boundaryFacets(shp);
                [x,y]=meshgrid(0:0.01:1,0:0.01:1);
                tri1=delaunay(x,y);
                tri=[];
                nodes=[];
                for i=1:size(BF,1)
                    plan1=[0 0 0;1 0 0;0 1 0]+2;
                    plan2=[P(BF(i,:),1),P(BF(i,:),2),P(BF(i,:),3)];
                    M=get_transform_plan(plan1,plan2);
                    nodes_plan=M*[x(:).'+2;y(:).'+2;ones(1,length(x(:)))*2;ones(1,length(x(:)))];
                    tri=[tri;tri1+length(nodes)];
                    nodes=[nodes;nodes_plan(1:3,:).'];
                end
                nodes=nodes*0.999;
                
                X0=sum([nodes(tri(:,1),1),nodes(tri(:,2),1),nodes(tri(:,3),1)],2)/3;
                Y0=sum([nodes(tri(:,1),2),nodes(tri(:,2),2),nodes(tri(:,3),2)],2)/3;
                Z0=sum([nodes(tri(:,1),3),nodes(tri(:,2),3),nodes(tri(:,3),3)],2)/3;
                
                in=shp.inShape(X0,Y0,Z0);
                tri=tri(in,:);
                X0=X0(in); Y0=Y0(in); Z0=Z0(in);
                w=wachspress([X0,Y0,Z0],obj);
                X=[nodes(tri(:,1),1),nodes(tri(:,2),1),nodes(tri(:,3),1)];
                Y=[nodes(tri(:,1),2),nodes(tri(:,2),2),nodes(tri(:,3),2)];
                Z=[nodes(tri(:,1),3),nodes(tri(:,2),3),nodes(tri(:,3),3)];
                c=permute(sum(couleur(1:n_couleurs,:).*(w.^5),1)+0.5,[3 2 1]);
                patch(X.',Y.',Z.',1:length(c),'CDataMapping','direct','Edgealpha',0);
                colormap(c); axis equal
                colorbar off; grid on;
                light('Position',[-50 -15 10])
                material([0.8 0.5 0.3])
                lighting flat
                view(45,30)
                
                for i=1:length(legend_mat)
                    text(v_center(i,1)*1.2,v_center(i,2)*1.2,v_center(i,3)*1.2,legend_mat(i))
                end
            end
            hold off;
            axis off
        end
        
        
        function plot_sub(obj,legende,k)
            color=["k","r","g","b","m","c","y"];
            if obj.dimension==0
                scatter3(obj.vertices(1),obj.vertices(2),obj.vertices(3),500,'filled',strcat(color(k),'p'));
                text(obj.vertices(1),obj.vertices(2),obj.vertices(3),legende)
                
            elseif obj.dimension==1
                
                plot3(obj.vertices(:,1),obj.vertices(:,2),obj.vertices(:,3),color(k),'linewidth',2)
                for i=1:length(legende)
                    text(obj.vertices(i,1),obj.vertices(i,2),obj.vertices(i,3),legende(i))
                end
            elseif obj.dimension==2
                obj.vertices(end+1,:)=obj.vertices(end,:)+1e-6;
                shp=alphaShape(obj.vertices,Inf);
                plot(shp,'facealpha',0.5,'facecolor',color(k),'edgecolor',color(k),'edgealpha',0.5);
                for i=1:length(legende)
                    text(obj.vertices(i,1),obj.vertices(i,2),obj.vertices(i,3),legende(i),'color',0.6*legend2color(legende(i))+0.5)
                end
            elseif obj.dimension==3
                shp=alphaShape(obj.vertices,Inf);
                plot(shp,'facealpha',0.5,'facecolor',color(k),'edgecolor',color(k),'edgealpha',0.5);
                
                for i=1:length(legende)
                    text(obj.vertices(i,1),obj.vertices(i,2),obj.vertices(i,3),legende(i),'color',0.6*legend2color(legende(i))+0.5)
                end
            end
        end
        
        
        function plot_normalFan(obj)
            alpha=0.2;
            L3D=1;
            L2D=3;
            h=figure();
            hold on
            listplot=obj.vertices(reshape(obj.edges.',[],1),:);
            
            if obj.dimension==2
                for i=1:length(obj.vertices) % point
                    p1=obj.vertices(i,:);
                    p2=p1+L2D*obj.normals(i,:);
                    p3=p1+L2D*obj.normals(mod(i+length(obj.vertices)-2,length(obj.vertices))+1,:);
                    p=[p1;p2;p3];
                    triangle=polyshape(p(:,1),p(:,2));
                    plot(triangle,'facecolor','r')
                end
                for i=1:length(obj.normal_fan.cones_edge) % edge
                    p1=[obj.vertices(obj.edges(i,1),1),obj.vertices(obj.edges(i,1),2)];
                    p4=[obj.vertices(obj.edges(i,2),1),obj.vertices(obj.edges(i,2),2)];
                    p2=p1+L2D*obj.normals(i,:);
                    p3=p4+L2D*obj.normals(i,:);
                    p=[p1;p2;p3;p4];
                    rectangle=polyshape(p(:,1),p(:,2));
                    plot(rectangle,'facecolor','g')
                end
                scatter(obj.vertices(:,1),obj.vertices(:,2),'k','filled');
                plot(reshape(listplot(:,1),2,[]),reshape(listplot(:,2),2,[]),'k');
            elseif obj.dimension==3
                for i=1:length(obj.normal_fan.cones_vertices) % point
                    fac_loc=obj.vertices2facets{i};
                    p=obj.vertices(i,:);
                    for j=1:length(fac_loc)
                        p=[p;p(1,:)+L3D*obj.normals(fac_loc(j),:)];
                    end
                    cone=alphaShape(p(:,1),p(:,2),p(:,3),Inf);
                    plot(cone,'facecolor','r','edgecolor','r','FaceAlpha',alpha,'EdgeAlpha',0)
                end
                for i=1:length(obj.normal_fan.cones_edges) % edge
                    fac_loc=obj.edges2facets(i,:);
                    p1=obj.vertices(obj.edges(i,1),:);
                    p2=obj.vertices(obj.edges(i,2),:);
                    p3=p1+obj.normals(fac_loc(1),:)*L3D;
                    p4=p2+obj.normals(fac_loc(1),:)*L3D;
                    p5=p1+obj.normals(fac_loc(2),:)*L3D;
                    p6=p2+obj.normals(fac_loc(2),:)*L3D;
                    p=[p1;p2;p3;p4;p5;p6];
                    toit=alphaShape(p(:,1),p(:,2),p(:,3),Inf);
                    plot(toit,'facecolor','g','edgecolor','g','FaceAlpha',alpha,'EdgeAlpha',0)
                end
                for i=1:length(obj.normal_fan.cones_facets) % facet
                    edges_facette=abs(obj.facets{i}(:));
                    nodes_facette=unique(obj.edges(edges_facette(:),:));
                    vloc_base=obj.vertices(nodes_facette,:);
                    vloc_extrude=vloc_base+obj.normals(i,:)*L3D;
                    vshp=[vloc_base;vloc_extrude];
                    prism=alphaShape(vshp(:,1),vshp(:,2),vshp(:,3),Inf);
                    plot(prism,'facecolor','b','edgecolor','b','FaceAlpha',alpha,'EdgeAlpha',0)
                end
                scatter3(obj.vertices(:,1),obj.vertices(:,2),obj.vertices(:,3),'k','filled');
                plot3(reshape(listplot(:,1),2,[]),reshape(listplot(:,2),2,[]),reshape(listplot(:,3),2,[]),'k');
                view(3)
            else
                error("incorrect dimension")
            end
            set(h,'renderer','Painters');
            axis equal
            axis off
            hold off
        end
        
    end
end

%% Other useful functions

function [rectangles,triangles,edges,vertices]=normalFan2D(v,epsilon)

if nargin<=2 || isempty(epsilon)
    epsilon=1e-7;
end

nv=length(v);
v=v-mean(v,1);
v=v*(1-epsilon);
edges=[(1:nv).',[(2:nv),1].'];
vertices=v;
normals=[0,1;-1 0]*(v(edges(:,2),:)-v(edges(:,1),:)).';
normals=((normals).')./vecnorm(normals.',2,2);

rectangles=cell(nv,1);

for i=1:nv % rectangle
    p1=[v(edges(i,1),1),v(edges(i,1),2)];
    p4=[v(edges(i,2),1),v(edges(i,2),2)];
    p2=p1+100*normals(i,:);
    p3=p4+100*normals(i,:);
    p=[p1;p2;p3;p4];
    rectangles{i}=polyshape(p(:,1),p(:,2));
end

triangles=cell(nv,1);

for i=1:nv-1 % triangles
    p1=v(i+1,:);
    p2=p1+100*normals(i,:);
    p3=p1+100*normals(i+1,:);
    p=[p1;p2;p3];
    triangles{i+1}=polyshape(p(:,1),p(:,2));
end
p1=v(1,:);
p2=p1+100*normals(nv,:);
p3=p1+100*normals(1,:);
p=[p1;p2;p3];
triangles{1}=polyshape(p(:,1),p(:,2));
end

function [prisms,roofs,cones,facets,edges,vertices]=normalFan3D(domain,epsilon)

if nargin<=2
    epsilon=1e-7;
end

v=domain.vertices;

v=v-mean(v,1);
v=v*(1-epsilon);

edges=domain.edges;
facets=domain.facets;
vertices=v;
normals=domain.normals;

L=100;

nv=length(v); ne=length(edges); nf=length(facets);

prisms=cell(nf,1);

for i=1:nf % prismes
    edges_facette=abs(facets{i}(:));
    nodes_facette=unique(edges(edges_facette(:),:));
    vloc_base=v(nodes_facette,:);
    vloc_extrude=vloc_base+normals(i,:)*L;
    vshp=[vloc_base;vloc_extrude];
    prisms{i}=alphaShape(vshp(:,1),vshp(:,2),vshp(:,3),Inf);
end

roofs=cell(nv,1);

for i=1:ne
    fac_loc=domain.edges2facets(i,:);
    p1=v(edges(i,1),:);
    p2=v(edges(i,2),:);
    p3=p1+normals(fac_loc(1),:)*L;
    p4=p2+normals(fac_loc(1),:)*L;
    p5=p1+normals(fac_loc(2),:)*L;
    p6=p2+normals(fac_loc(2),:)*L;
    p=[p1;p2;p3;p4;p5;p6];
    roofs{i}=alphaShape(p(:,1),p(:,2),p(:,3),Inf);
end

cones=cell(nv,1);

for i=1:nv % cones
    fac_loc=domain.vertices2facets{i};
    p=v(i,:);
    for j=1:length(fac_loc)
        p=[p;p(1,:)+L*normals(fac_loc(j),:)];
    end
    cones{i}=alphaShape(p(:,1),p(:,2),p(:,3),Inf);
end

end

function out1= get_transform_plan(plane1,plane2)

x1=plane1(1,1);
x2=plane1(2,1);
x3=plane1(3,1);

y1=plane1(1,2);
y2=plane1(2,2);
y3=plane1(3,2);

z1=plane1(1,3);
z2=plane1(2,3);
z3=plane1(3,3);

a1=plane2(1,1);
a2=plane2(2,1);
a3=plane2(3,1);

b1=plane2(1,2);
b2=plane2(2,2);
b3=plane2(3,2);

c1=plane2(1,3);
c2=plane2(2,3);
c3=plane2(3,3);

t2 = x1.*y2.*z3;
t3 = x1.*y3.*z2;
t4 = x2.*y1.*z3;
t5 = x2.*y3.*z1;
t6 = x3.*y1.*z2;
t7 = x3.*y2.*z1;
t8 = -t3;
t9 = -t4;
t10 = -t7;
t11 = t2+t5+t6+t8+t9+t10;
t12 = 1.0./t11;
mt1 = [t12.*(a1.*y2.*z3-a1.*y3.*z2-a2.*y1.*z3+a2.*y3.*z1+a3.*y1.*z2-a3.*y2.*z1),t12.*(b1.*y2.*z3-b1.*y3.*z2-b2.*y1.*z3+b2.*y3.*z1+b3.*y1.*z2-b3.*y2.*z1),t12.*(c1.*y2.*z3-c1.*y3.*z2-c2.*y1.*z3+c2.*y3.*z1+c3.*y1.*z2-c3.*y2.*z1),t12.*(y1.*z2-y2.*z1-y1.*z3+y3.*z1+y2.*z3-y3.*z2),-t12.*(a1.*x2.*z3-a1.*x3.*z2-a2.*x1.*z3+a2.*x3.*z1+a3.*x1.*z2-a3.*x2.*z1),-t12.*(b1.*x2.*z3-b1.*x3.*z2-b2.*x1.*z3+b2.*x3.*z1+b3.*x1.*z2-b3.*x2.*z1)];
mt2 = [-t12.*(c1.*x2.*z3-c1.*x3.*z2-c2.*x1.*z3+c2.*x3.*z1+c3.*x1.*z2-c3.*x2.*z1),-t12.*(x1.*z2-x2.*z1-x1.*z3+x3.*z1+x2.*z3-x3.*z2),t12.*(a1.*x2.*y3-a1.*x3.*y2-a2.*x1.*y3+a2.*x3.*y1+a3.*x1.*y2-a3.*x2.*y1),t12.*(b1.*x2.*y3-b1.*x3.*y2-b2.*x1.*y3+b2.*x3.*y1+b3.*x1.*y2-b3.*x2.*y1),t12.*(c1.*x2.*y3-c1.*x3.*y2-c2.*x1.*y3+c2.*x3.*y1+c3.*x1.*y2-c3.*x2.*y1),t12.*(x1.*y2-x2.*y1-x1.*y3+x3.*y1+x2.*y3-x3.*y2),0.0,0.0,0.0,0.0];
out1 = reshape([mt1,mt2],4,4);
end




function rho=projection2D(rho,triangles,rectangles,vertices,edges,epsilon)

if nargin <=5 || isempty(epsilon)
    epsilon=1e-7;
end

% 1) projection onto vertices (triangles)
N=length(triangles);
for i=1:N
    s=triangles{i}.Vertices;
    in=inpolygon(rho(:,1),rho(:,2),s(:,1),s(:,2));
    if sum(in)>0
        rho(in,1)=vertices(i,1)*(1-epsilon);
        rho(in,2)=vertices(i,2)*(1-epsilon);
    end
end

% 2) projection onto edges (rectangles)

for i=1:N
    s=rectangles{i}.Vertices;
    in=inpolygon(rho(:,1),rho(:,2),s(:,1),s(:,2));
    if sum(in)>0
        projpoints=proj([vertices(edges(i,:).',1),vertices(edges(i,:).',2)]*(1-epsilon), rho(in,:));
        rho(in,:)=shiftdim(projpoints,2);
    end
end


end

function [ProjPoint] = proj(vector, q)

p0 = vector(1,:);
p1 = vector(2,:);
q=shiftdim(q,-2);

a = [-q(1,1,:,1).*(p1(1)-p0(1)) - q(1,1,:,2).*(p1(2)-p0(2)); ...
    repmat(-p0(2).*(p1(1)-p0(1)) + p0(1).*(p1(2)-p0(2)),[ 1,1,size(q,3)])];
b = [p1(1) - p0(1), p1(2) - p0(2);...
    p0(2) - p1(2), p1(1) - p0(1)];

detb=b(1,1).*b(2,2)-b(1,2).*b(2,1);
invb=[b(1,1), b(2,1); b(1,2), b(2,2)]./detb;

ProjPoint = -mult(invb,a);
end

function rho=projection3D(rho,un,prisms,roofs,cones,facets,edges,vertices,epsilon)

if nargin<=8 || isempty(epsilon)
    epsilon=1e-7;
end
% 1) projection onto vertices
N=length(cones);
vertices=vertices-mean(vertices,1);
for i=1:N
    in=inShape(cones{i},rho);
    if sum(in)>0
        rho(in,1)=(1-epsilon)*vertices(i,1);
        rho(in,2)=(1-epsilon)*vertices(i,2);
        rho(in,3)=(1-epsilon)*vertices(i,3);
    end
end

% 2) projection onto edges

N=length(roofs);
for i=1:N
    in=inShape(roofs{i},rho);
    if sum(in)>0
        p=vertices(edges(i,1),:);
        vec_dir=vertices(edges(i,2),:)-vertices(edges(i,1),:);
        projpoints=proj_edge((1-epsilon)*vec_dir,(1-epsilon)*p,rho(in,:));
        rho(in,:)=shiftdim(projpoints,2);
    end
end

% 3) projection onto facets

N=length(prisms);
for i=1:N
    in=inShape(prisms{i},rho);
    if sum(in)>0
        p=vertices(edges(abs(facets{i}(1)),1),:);
        n=un(i,:);
        projpoints=proj_plan(n,(1-epsilon)*p,rho(in,:));
        rho(in,:)=shiftdim(projpoints,2);
    end
end

end

function [ProjPoint] = proj_edge(vec_dir,p,q)

q=permute(q,[3,2,1]);
vec_dir=vec_dir(:)/norm(vec_dir,2);
ProjPoint=(p+mult(q-p,vec_dir).*vec_dir.');
end

function [ProjPoint] = proj_plan(n,p, q)
% project onto the plane
q=permute(q,[3,2,1]);
ProjPoint = (q - mult(q - p, n(:)).*n);
end

function plotcolor2D(p,elts,data)
elt1=elts(:,1);
elt2=elts(:,2);
elt3=elts(:,3);
x=[p(elt1,1), p(elt2,1), p(elt3,1)];
y=[p(elt1,2), p(elt2,2), p(elt3,2)];
z = zeros(size(x));
patch(x.',y.',z.',data(:).','Edgecolor','none');
end

function colors=legend2color(str_legend,inv)

% associate a color to a string

if nargin<=1 || isempty(inv)
    inv=0;
end

dict_colors.fer=[0 0 0]-0.5;
dict_colors.iron=[0 0 0]-0.5;
dict_colors.fesi=[0 0 0]-0.5;
dict_colors.fesi1=[0 0 0]-0.5;
dict_colors.fesi2=[0 0 0]-0.5;
dict_colors.feco=[0 0 0]-0.5;
dict_colors.feni=[0 0 0]-0.5;
dict_colors.smc=[0 0 0]-0.5;

dict_colors.air=[1 1 1]-0.5;

dict_colors.ap=[1 0 0]-0.5;
dict_colors.am=[1 0 1]-0.5;
dict_colors.bp=[0 1 0]-0.5;
dict_colors.bm=[1 1 0]-0.5;
dict_colors.cp=[0 0 1]-0.5;
dict_colors.cm=[0 1 1]-0.5;


dict_colors.("v1")=[0 0 0]-0.5;
dict_colors.("v2")=[1 0 0]-0.5;
dict_colors.("v3")=[1 0 1]-0.5;
dict_colors.("v4")=[0 1 0]-0.5;
dict_colors.("v5")=[1 1 0]-0.5;
dict_colors.("v6")=[0 0 1]-0.5;
dict_colors.("v7")=[0 1 1]-0.5;
dict_colors.("v8")=[1 1 1]-0.5;
dict_colors.("v9")=[0.8 0.2 0.2]-0.5;
dict_colors.("v10")=[0.8 0.2 0.8]-0.5;
dict_colors.("v11")=[0.2 0.8 0.2]-0.5;
dict_colors.("v12")=[0.8 0.8 0.2]-0.5;
dict_colors.("v13")=[0.2 0.2 0.8]-0.5;
dict_colors.("v14")=[0.2 0.8 0.8]-0.5;
dict_colors.("v15")=[0.6 0.1 0.1]-0.5;
dict_colors.("v16")=[0.1 0.1 0.6]-0.5;
dict_colors.("v17")=[0.1 0.6 0.1]-0.5;
dict_colors.("v18")=[0.6 0.6 0.1]-0.5;
dict_colors.("v19")=[0.1 0.1 0.6]-0.5;
dict_colors.("v20")=[0.1 0.1 0.6]-0.5;

dict_colors.dp=[0.8 0.2 0.2]-0.5;
dict_colors.dm=[0.8 0.2 0.8]-0.5;
dict_colors.ep=[0.2 0.8 0.2]-0.5;
dict_colors.em=[0.8 0.8 0.2]-0.5;
dict_colors.fp=[0.2 0.2 0.8]-0.5;
dict_colors.fm=[0.2 0.8 0.8]-0.5;
dict_colors.gp=[0.6 0.1 0.1]-0.5;
dict_colors.gm=[0.1 0.1 0.6]-0.5;
dict_colors.hp=[0.1 0.6 0.1]-0.5;
dict_colors.hm=[0.6 0.6 0.1]-0.5;
dict_colors.ip=[0.1 0.1 0.6]-0.5;
dict_colors.im=[0.1 0.1 0.6]-0.5;

dict_colors.mag1=[1 0 0]-0.5;
dict_colors.mag2=[1 0 1]-0.5;
dict_colors.mag3=[0 1 0]-0.5;
dict_colors.mag4=[1 1 0]-0.5;
dict_colors.mag5=[0 0 1]-0.5;
dict_colors.mag6=[0 1 1]-0.5;
dict_colors.mag7=[1 0 0]-0.5;
dict_colors.mag8=[1 0 1]-0.5;
dict_colors.mag9=[0 1 0]-0.5;
dict_colors.mag10=[1 1 0]-0.5;
dict_colors.mag11=[0 0 1]-0.5;
dict_colors.mag12=[0 1 1]-0.5;

dict_colors.magnet1=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet2=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet3=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet4=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet5=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet6=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet7=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet8=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet9=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet10=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet11=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet12=[0.3 0.3  0.3 ]-0.5;


dict_colors.mag_h=[1 0 0]-0.5;
dict_colors.mag_b=[1 0 1]-0.5;
dict_colors.mag_l=[0 1 0]-0.5;
dict_colors.mag_r=[1 1 0]-0.5;

dict_colors.mag_n=[1 0 0]-0.5;
dict_colors.mag_s=[1 0 1]-0.5;
dict_colors.mag_e=[0 1 0]-0.5;
dict_colors.mag_w=[1 1 0]-0.5;

dict_colors.mag_up=[1 0 0]-0.5;
dict_colors.mag_down=[1 0 1]-0.5;
dict_colors.mag_right=[0 1 0]-0.5;
dict_colors.mag_left=[1 1 0]-0.5;


dict_colors.fer1=[0 0 0]-0.5;
dict_colors.fer2=[0 0 1]-0.5;
dict_colors.fer3=[0 1 0]-0.5;
dict_colors.fer4=[0 1 1]-0.5;
dict_colors.fer5=[1 0 0]-0.5;
dict_colors.fer6=[1 0 1]-0.5;
dict_colors.fer7=[1 1 0]-0.5;

dict_colors.conducteur=[1 0.5 0]-0.5;
dict_colors.aimant=[1 0 0]-0.5;
dict_colors.conducteurs=[1 0.5 0]-0.5;
dict_colors.aimants=[1 0 0]-0.5;
dict_colors.mag=[1 0 0]-0.5;
dict_colors.fers=[1 0 0]-0.5;

dict_colors.inductor=[0.8 0 0.2]-0.5;
dict_colors.inducteur=[0.8 0 0.2]-0.5;
dict_colors.indp=[0.8 0 0.2]-0.5;
dict_colors.indm=[0.2 0 0.8]-0.5;


if nargin==0 || isempty(str_legend)
    colors=convertCharsToStrings(fieldnames(dict_colors));
else
    str_legend(lower(str_legend)=="a+")="ap";
    str_legend(lower(str_legend)=="a-")="am";
    str_legend(lower(str_legend)=="b+")="bp";
    str_legend(lower(str_legend)=="b-")="bm";
    str_legend(lower(str_legend)=="c+")="cp";
    str_legend(lower(str_legend)=="c-")="cm";
    str_legend(lower(str_legend)=="d+")="dp";
    str_legend(lower(str_legend)=="d-")="dm";
    str_legend(lower(str_legend)=="e+")="ep";
    str_legend(lower(str_legend)=="e-")="em";
    str_legend(lower(str_legend)=="f+")="fp";
    str_legend(lower(str_legend)=="f-")="fm";
    str_legend(lower(str_legend)=="g+")="gp";
    str_legend(lower(str_legend)=="g-")="gm";
    str_legend(lower(str_legend)=="h+")="hp";
    str_legend(lower(str_legend)=="h-")="hm";
    str_legend(lower(str_legend)=="i+")="ip";
    str_legend(lower(str_legend)=="i-")="im";
    
    str_legend(lower(str_legend)=="ind-")="indm";
    str_legend(lower(str_legend)=="ind+")="indp";
    
    legend2=str_legend;
    if inv
        legend2(lower(str_legend)=="ap")="am";
        legend2(lower(str_legend)=="am")="ap";
        legend2(lower(str_legend)=="bm")="bp";
        legend2(lower(str_legend)=="bp")="bm";
        legend2(lower(str_legend)=="cm")="cp";
        legend2(lower(str_legend)=="cp")="cm";
        legend2(lower(str_legend)=="dp")="dm";
        legend2(lower(str_legend)=="dm")="dp";
        legend2(lower(str_legend)=="em")="ep";
        legend2(lower(str_legend)=="ep")="em";
        legend2(lower(str_legend)=="fm")="fp";
        legend2(lower(str_legend)=="fp")="fm";
        legend2(lower(str_legend)=="gp")="gm";
        legend2(lower(str_legend)=="gm")="gp";
        legend2(lower(str_legend)=="hm")="hp";
        legend2(lower(str_legend)=="hp")="hm";
        legend2(lower(str_legend)=="im")="ip";
        legend2(lower(str_legend)=="ip")="im";
        
        legend2(lower(str_legend)=="mag_up")="mag_down";
        legend2(lower(str_legend)=="mag_down")="mag_up";
        legend2(lower(str_legend)=="mag_n")="mag_s";
        legend2(lower(str_legend)=="mag_s")="mag_n";
        legend2(lower(str_legend)=="mag_e")="mag_w";
        legend2(lower(str_legend)=="mag_w")="mag_e";
        legend2(lower(str_legend)=="mag_right")="mag_left";
        legend2(lower(str_legend)=="mag_left")="mag_right";
        legend2(lower(str_legend)=="mag_h")="mag_b";
        legend2(lower(str_legend)=="mag_b")="mag_h";
        legend2(lower(str_legend)=="mag_l")="mag_r";
        legend2(lower(str_legend)=="mag_r")="mag_l";
        
        legend2(lower(str_legend)=="indp")="indm";
        legend2(lower(str_legend)=="indm")="indp";
    end
    
    
    for i=1:length(legend2)
        try
            colors(i,:)=dict_colors.(lower(legend2(i)));
        catch
            colors(i,:)=[0,0,0];
        end
    end
    
end


end