
    clear;
    %%%%close all;
    clc;
    
    kernels = double(importdata('allkernels.mat')); 
    shells = double(importdata('allshells.mat')); 
    
    [x,y] = meshgrid( ...
        linspace(0, 1, 100), ...
        linspace(0, 1, 100) ...
    );

    [rows, cols] = size(x);

    doe = DoEstimator();
    
    
    colormap hsv;
    idx = 1;
    for deg = 0:0.05:0.45
        n = size(kernels, 1);
        disp(deg);
        contamination = shells(randi(size(shells, 1), floor(n * deg), 1), :);
        dataset = [kernels; contamination];
        [dotrain, dotest] = doe.trainModel(dataset(:, [1 2]), [x(:), y(:)]);        
        %%subplot(3,4, idx);
        figure; 
        contour(x,y, reshape(dotest, [cols, rows]), 100); shading flat; hold on;        
        plot3(dataset(:, 1), dataset(:, 2), dotrain, '.'); colorbar; title([mat2str(100 * deg) '% contamination added']);
        idx = idx+1;
    end
    
    
    
    
    
    
    
    
    