function stop = outfunc( t, y, state , varargin )
stop=false;
persistent step
persistent t0
persistent tt0
spacing = 1e4;

tt=toc;

switch state
    case []
        step = step + numel(t); 

        if mod(step,spacing)==0
            
            dt = (t-t0)/spacing;
            dtt = (tt-tt0);
            disp(['Step ' num2str(step) ', wall ' num2str(tt) ', sim ' num2str(t(1))...
            	  ', s/w ' num2str(dt(1)/dtt) ', avg step ' num2str(dt(1))]);
            t0 = t;
            tt0 = tt;
        end

    case 'init'
        step = 0;
        t0=0;
        tt0=0;
        disp(' ');
end
end
