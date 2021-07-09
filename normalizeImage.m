function nM = normalizeImage(img)

nM = (img - min(img(:)))/(max(img(:)) - min(img(:)));