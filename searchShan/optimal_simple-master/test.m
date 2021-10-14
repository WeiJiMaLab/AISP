x = -5; sigma = 5; sigma_s = 9.0593;
sTMat = linspace(0,100,10000);

anlT = sqrt(2*pi/(1/sigma^2+1/sigma_s^2))*(1/2+1/2*erf(x/sigma^2/sqrt(2*(1/sigma^2+1/sigma_s^2))))*exp(-x^2/2/(sigma^2+sigma_s^2))

numT = sum(exp(-(x-sTMat).^2/2/sigma^2).*exp(-sTMat.^2/2/sigma_s^2))*range(sTMat)/length(sTMat)


xD = [-5,-5,8]; N = 4;
sDMat = linspace(-100,100,10000);

weight = (sum(xD)^2)/(sigma^4)/((N-1)/sigma^2+1/sigma_s^2) - sum(xD.^2)/sigma^2;
anlD = sqrt(2*pi/((N-1)/sigma^2+1/sigma_s^2))*exp(weight/2)

temp = 1;
for ii = 1:N-1
    temp = temp.*exp(-(xD(ii)-sDMat).^2/2/sigma^2);
end

xDMat = repmat(xD',1,length(sDMat));
temp2 = squeeze(prod(exp(-bsxfun(@minus,xDMat,sDMat).^2/2/sigma^2)));

numD = sum(temp.*exp(-sDMat.^2/2/sigma_s^2))*range(sDMat)/length(sDMat)
numD2 = sum(temp2.*exp(-sDMat.^2/2/sigma_s^2))*range(sDMat)/length(sDMat)


weight1 = (sum(xD)^2)/(sigma^4)/((N-1)/sigma^2+1/sigma_s^2) - sum(xD.^2)/sigma^2
weight2 = mean(xD)^2/(sigma_s^2+sigma^2/(N-1)) + var(xD,1)*(N-1)/sigma^2


weight3 = (sum(xD)^2)/(sigma^4)/((N-1)/sigma^2+1/sigma_s^2) - sum(xD.^2)/sigma^2
weight4 = sum(xD)^2/(N-1)/sigma^2 - sum(xD)^2/sigma^2/sigma_s^2/((N-1)/sigma^2+1/sigma_s^2)/(N-1)- sum(xD.^2)/sigma^2




% x = [-5,-5,-5,4];
% term = zeros(1,N);weight = zeros(1,N); idx = 1:N;
% for ii = 1:N
%     xT = x(idx==ii); xD = x(idx~=ii);
%     weight(ii) = xT^2/(sigma^2+sigma_s^2) + mean(xD)^2/(sigma_s^2+sigma^2/(N-1)) + var(xD)*(N-1)/sigma^2;
%     term(ii) = (1/2+1/2*erf(xT/sigma^2/sqrt(2*(1/sigma^2+1/sigma_s^2))));
%     anal = sum(term(ii).*exp(-weight(ii)/2))*sqrt(2*pi/((N-1)/sigma^2+1/sigma_s^2))*sqrt(2*pi/(1/sigma^2+1/sigma_s^2))
%     term1 = sum(exp(-(xT-sTMat).^2/2/sigma^2).*exp(-sTMat.^2/2/sigma_s^2))*range(sTMat)/length(sTMat);
%     temp = 1;
%     for jj = 1:N-1
%         temp = temp.*exp(-(xD(jj)-sDMat).^2/2/sigma^2);
%     end
%     term2 = sum(temp.*exp(-sDMat.^2/2/sigma_s^2))*range(sDMat)/length(sDMat);
%     num = term1*term2
% end
% 
