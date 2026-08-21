% function to find the center of beads in the nine panels with 2d gaussian
% fitting.
% written by Bassam Hajj
% 31/12/2013

%input: Max intensity projection image
%       starting location (x,y)
%       aa : window's size
%output: coordinates


function coordinates=gaussian_fitting(MIP_image,start_location,aa)



[x y]=meshgrid(0:2*aa,0:2*aa);
fitparams=zeros(9,5);
        
        
        
        
            %tt=[start_location(1)-aa:start_location(1)+aa;start_location(2)-aa:start_location(2)+aa];
            patch2=MIP_image(start_location(2)-aa:start_location(2)+aa,start_location(1)-aa:start_location(1)+aa);

            [r c] =find(patch2==max(patch2(:))); %finds the maximum value in the matrix
            %recentring arround this value
            start_location(1)= max(c)+start_location(1)-aa;
            start_location(2)= max(r)+start_location(2)-aa;

            patch2=MIP_image(start_location(2)-aa:start_location(2)+aa,start_location(1)-aa:start_location(1)+aa);
            meanbackground=mean([mean(patch2(1,:)) mean(patch2(2*aa+1,:)) mean(patch2(:,1)) mean(patch2(:,2*aa+1))]);
           % patch2=patch2-meanbackground;% remooving background
            %define the Gaussian 2D function
            myfun = @(A) A(1)*exp(-(((x(:)-A(2)).^2/(2*A(3)^2))+((y(:)-A(4)).^2/(2*A(5)^2))))-patch2(:);

            [r c] =find(patch2==max(patch2(:)));
            % Initial guess parameters
            r=r(1);
            c=c(1);
           
            A0=[patch2(r,c);aa;2;aa;2];
            %Add lots of debugging info
            %opts = optimset('Display','Iter');
            options = optimset('Algorithm',{'levenberg-marquardt',0.005},'display','off');
            fitparams = (lsqnonlin(myfun,A0,[],[],options))';      
        
            
            
            coordinates(1)=fitparams(2)+start_location(1)-aa;
            coordinates(2)=fitparams(4)+start_location(2)-aa;
           % figure(10),imagesc(patch2),hold on,plot(fitparams(4),fitparams(2),'rx')

end