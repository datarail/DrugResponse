function savegcf(filename,opengl)
% savegcf(filename)

if exist('opengl', 'var') && opengl 
    set(gcf,'Renderer','opengl')
else
    set(gcf,'Renderer','painters')
end
if nargin==0 || isempty(filename)
    saveas(gcf,get(gcf,'FileName'))
else
    saveas(gcf,filename)
end
