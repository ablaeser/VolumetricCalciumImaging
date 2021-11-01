%function to transfer
fun = 'C:\Users\LehtinenLab\Dropbox\AndermannLab\users\Fred\pipe-master\Fbs registration pipeline\SmallAreaRegistration.m';
%destination directory
dest = 'C:\Users\LehtinenLab\Dropbox\Jin\MIA project Figures\Cui2020-code';
%get all dependencies
[fList,pList] = matlab.codetools.requiredFilesAndProducts(fun);
%copy each dependency to the directory
for i = 1:numel(fList)
    copyfile(fList{i},dest);
end