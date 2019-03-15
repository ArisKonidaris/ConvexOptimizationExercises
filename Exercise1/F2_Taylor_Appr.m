function [FirstTaylor, SecondTaylor] = F2_Taylor_Appr(Num1,Num2,X,Y)
  FirstTaylor = 1.-((Num1+Num2)^2.+X.+Y)./(1+Num1+Num2)^2;
  SecondTaylor = FirstTaylor.+((X.-Num1.+Y.-Num2).^2)...
                 ./(1+Num1+Num2)^3;
  return
endfuncion