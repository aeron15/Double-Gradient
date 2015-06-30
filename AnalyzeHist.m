function [perc_const,perc_adaptive,perc_area] = AnalyzeHist(y,x,th_const,off_peak,off_dist_y,off_dist_x);

p = y./sum(y);

if ~isempty(find(x>=th_const))
    ind= find(x>=th_const);
    perc_const = sum(y(ind:end))/sum(y(1:ind-1));
    perc_adaptive=perc_const;
else
    perc_const = 0;
    perc_adaptive=0;
end

% set up a common x 

max_x = max([max(x),max(off_dist_x)]);
min_x = min([min(x),min(off_dist_x)]);
xx = linspace(min(x),max(x),150);

yy = spline(x,y,xx);yy(yy<0)=0;
yy_off = spline(off_dist_x,off_dist_y,xx);yy_off(yy_off<0)=0;

pp = yy/sum(yy);
p_off =yy_off/sum(yy_off);
diff_p = pp-p_off;diff_p(diff_p<0)=0;
perc_area = sum(diff_p);