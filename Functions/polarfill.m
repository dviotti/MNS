function polarfill(ax_polar,theta,rlow,rhigh,color,alpha)
    ax_cart = axes();
    ax_cart.Position = ax_polar.Position;
    [xl,yl] = pol2cart(theta,rlow);
    [xh,yh] = pol2cart(fliplr(theta),fliplr(rhigh));
    fill([xl,xh],[yl,yh],color,'FaceAlpha',alpha,'EdgeAlpha',0);
    xlim(ax_cart,[-max(get(ax_polar,'RLim')),max(get(ax_polar,'RLim'))]); 
    ylim(ax_cart,[-max(get(ax_polar,'RLim')),max(get(ax_polar,'RLim'))]);
    axis square; set(ax_cart,'visible','off');
end