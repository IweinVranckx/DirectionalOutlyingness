classdef DoEstimator
    
    properties
        
    end
    
    methods (Access = private)
        function r=rho(~, x)
            c = 2.1;
            r=(x / c).^2;
            r(r>1)=1;
            r = 1.54^2*r;
        end
        
        function result = DO(~, projectedSamples, mu, Sa, Sb)
           assert(isvector(projectedSamples), 'samples should be univariate');
           
           result = zeros(length(projectedSamples), 1);
           mask = projectedSamples>=mu;
           result(mask) = (projectedSamples(mask) - mu) ./ Sa;
           result(~mask) = (projectedSamples(~mask) - mu) ./ Sb;
        end
        
    end
    
    methods     
        function this = DoEstimator()        
        end
        
        %%%%    https://www.math.washington.edu/~king/coursedir/m445w04/notes/vector/equations.html
        %%%%    https://www.math.washington.edu/~king/coursedir/m445w04/notes/vector/normals-planes.html
        function [doTraining, doTest] = trainModel(this, trainSamples, testSamples)
            
            assert(nargin==3, '3 arguments needed');
            
            [n, p] = size(trainSamples);
            nbrOfdirections = 250 * p;
            
            h = floor((n+1)/2);            
            b = ones(p, 1);
            
            doTraining = zeros(size(trainSamples, 1), 1);
            doTest = zeros(size(testSamples, 1), 1);
            
            for dindex=1:nbrOfdirections
                
                %%%%    Train
                
                subset = trainSamples(randi([1, n], [p, 1]), :);
                normalVector = linsolve(subset, b);
                normalVector = normalVector / sqrt(sum(normalVector.^2));
                prj = trainSamples * normalVector;
                
                Y=sort(prj);        
                if rem(length(prj), 2)    % It's odd
                    Ya = Y(h:end);
                    Yb = Y(1:h);
                else
                    Ya = Y(h+1:end);
                    Yb = Y(1:h);
                end

                mu = median(prj);
                Za = Ya  - mu;
                Zb = mu - Yb;
                Soa = 1.4826 * median(Za);
                Sob = 1.4826 * median(Zb);
                %%% one step M-estimates are:    
                Sa = Soa * sqrt((2/length(Za)) * sum(this.rho(Za/Soa)));
                Sb = Sob * sqrt((2/length(Zb)) * sum(this.rho(Zb/Sob)));
                
                doTraining = max(doTraining, this.DO(prj, mu, Sa, Sb));
                
                %%% Het schijnt niet mogelijk te zijn om de beste
                %%% projectievector zomaar te bepalen uit 'n' voorstellen.
                
                %%%%    Test phase                
                projTest = normalVector' * testSamples';
                doTest = max(doTest, this.DO(projTest, mu, Sa, Sb));
            end
        end
    end
    
end

