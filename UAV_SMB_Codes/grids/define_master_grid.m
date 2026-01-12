function [Xq, Yq] = define_master_grid(R, target_res)
xmin = R.XWorldLimits(1);
xmax = R.XWorldLimits(2);
ymin = R.YWorldLimits(1);
ymax = R.YWorldLimits(2);

x_vec = xmin:target_res:xmax;
y_vec = ymin:target_res:ymax;

[Xq, Yq] = meshgrid(x_vec, y_vec);
end
