%delete('if-neonate-spontaneous.zip');
for i = 1:length(zipfilecontents)
    delete(zipfilecontents{i});
end
