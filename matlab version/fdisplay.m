function [] = fdisplay(X,Y,A)
% h=surf(X,Y,A);
h = contourf(X,Y,A);
view(0,90);
% set(h,'LineStyle','none');
end