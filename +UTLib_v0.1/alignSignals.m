function tsense = alignSignals(texcite, yexcite, tsense)

[~, I] = max(yexcite);
tshift = texcite(I);
tsense = tsense+tshift;

end

