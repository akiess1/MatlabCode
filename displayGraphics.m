        

function displayGraphics(xmin,xmax,Nx,ymin,ymax,Ny,bacteria,alive,c1,c2,clrscl,...
     t,figHand)

    dx = (xmax - xmin)/Nx;
    dy = (ymax - ymin)/Ny;
    x = xmin:dx:xmax; % including ghost nodes x(1) and x(Nx+3)
    y = ymin:dy:ymax; % including ghost nodes y(1) and y(Ny+3)
    
    if alive == 1 
        Nb = length(bacteria);
        colors = zeros(Nb,3);
        xpos = zeros(1,Nb);
        ypos = zeros(1,Nb);
        for bi = 1:Nb
            xpos(bi) = bacteria(bi).position(1);
            ypos(bi) = bacteria(bi).position(2);
            colors(bi,:) = (bacteria(bi).secretion(1)*[0 0.01 0])+[0.1 0.2 0];
        end
    end
    
    Clrs = zeros(Ny+1,Nx+1,3);
    for xi = 2:Nx+2
        for yi = 2:Ny+2
            Clrs(yi-1,xi-1,:) = [clrscl 0 0]*c1(xi,yi) ...
                + [0 0 clrscl]*c2(xi,yi);
        end
    end        
    
    figure(figHand);
    image(x,y,Clrs);
    set(gca,'YDir','normal');
    shading interp;
    hold on;
    if alive ==1
        scatter(xpos,ypos,8,colors,'filled');
    end
    axis([xmin xmax ymin ymax]);

    title(sprintf('time = %1.8f',t)); 
    pause(0.0001);
    hold off;
    
        %{
        if saveVid == 1
            if mod(ti,saveRate) == 1;
               F = getframe(gcf);
               writeVideo(vid,F);
            end
        end
        %}
 
end

