% Generate channNum channel instantiations of mmWave channel 
% By Le Liang, UVic, July 6, 2014

function [genH, genAlpha, genAt, genAr] = channelSet(Nt, Nr, Ncls, Nray, channNum, DirTX, DirRX, large_fa_loss)

% display('Generating channel...')

genH = zeros(prod(Nr), prod(Nt), channNum);
genAlpha = zeros(Ncls*Nray, channNum);
genAt = zeros(prod(Nt), Ncls*Nray, channNum);
genAr = zeros(prod(Nr), Ncls*Nray, channNum);

for ichann = 1 : channNum
    [Ht, alphat, Att, Art] = GenChannel(Nt, Nr, Ncls, Nray, DirTX, DirRX);
    % add large fading loss ==> sqrt(large_fa_loss) ~ 1
    loss = (rand() * (1 - sqrt(large_fa_loss)) + sqrt(large_fa_loss));
    genH(:, :, ichann) = Ht * loss;  
    genAlpha(:, ichann) = alphat * loss;
    genAt(:, :, ichann) = Att;
    genAr(:, :, ichann) = Art;
end

pathname = fileparts('ChannelData\');% specify the pathname

% outfile = fullfile(pathname, sprintf('DIRchannel-Nt%d-Nr%d-Ncls%d-Nray%d-channNum%d',prod(Nt),prod(Nr),Ncls,Nray,channNum));
% save(outfile,'genH','genAlpha','genAt','genAr') % save channel data

outfile = fullfile(pathname, sprintf('OMNIchannel-Nt%d-Nr%d-Ncls%d-Nray%d-channNum%d',prod(Nt),prod(Nr),Ncls,Nray,channNum));
% save(outfile,'genH','genAlpha','genAt','genAr'); % save channel data

% display('Done!')
