function DV = SEM_DV(img)
dip = dipesti2D(hilbert(img));
DV = dipSemb2D03(dip,img);
DV = hist(DV(:),100);
end

