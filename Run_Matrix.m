%Run Matrix

function Run_Matrix()
clear all; clc; close all;

    load RunMatrix.mat

      for i = 1:length(RunMatrix.RunOrder)
        tic
        [fval,exitflag,loop,maxStress,xnew,volfrac] = ...
            mainMultipleLoads(RunMatrix.volfrac(i),RunMatrix.penalty(i),RunMatrix.r_min(i));
        i
        time = toc;
        RunMatrix.fval(i,1) = fval;
        RunMatrix.volfrac_final(i,1) = volfrac;
        RunMatrix.exitflag(i,1) = exitflag;
        RunMatrix.loop(i,1) = loop;
        RunMatrix.maxStress(i,1) = maxStress;
        RunMatrix.time(i,1) = time;
        
        x(:,i) = xnew;     
      end

save('multipleLoadsData')

end