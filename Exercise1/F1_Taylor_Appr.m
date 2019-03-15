function [FirstTaylor, SecondTaylor] = F1_Taylor_Appr(Num,Dom)
  FirstTaylor = (1-Dom+2*Num)./(1+Num)^2;
  SecondTaylor = (1+3*Num^2+3*Num-3*Num.*Dom.-Dom.+Dom.^2)...
                 ./(1+Num)^3;
  return
endfuncion