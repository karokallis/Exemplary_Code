function mb0 = estimateCalibratedmb0(prediction,f)
%% Created by Karoline Kallis 07/2023
% Method to synth. DWI b=0
% calculate mb0 scalar for each patient
% estimate correction Factor between different Sequences

load('./data/params.mat');
bvals = prediction.data.bvalues;
rsi_predicted = zeros([size(prediction.data.outmaps{1}) length(prediction.data.outmaps)]);
ProstateSegmentation = zeros([size(prediction.data.outseg{1}) length(prediction.data.outmaps)]);
calcDWI = zeros(size(rsi_predicted,1),size(rsi_predicted,2),size(rsi_predicted,3),numel(bvals),size(rsi_predicted,5));
mb0 = zeros(length(prediction.data.outmaps),1);

%% Calcualte calibrated mb0
for pati=1:length(prediction.data.outmaps)
    if any(size(prediction.data.outseg{pati})~=size(ProstateSegmentation(:,:,:,pati)))
       ProstateSegmentation(:,:,:,pati) = imresize3(prediction.data.outseg{pati},size(ProstateSegmentation(:,:,:,pati))); 
    else
        ProstateSegmentation(:,:,:,pati) = prediction.data.outseg{pati};
    end
    for bvali=1:numel(bvals)
        for cvali=1:size(rsi_predicted,4)
            %              rsi_predicted(:,:,:,cvali,pati) = prediction.data.outmaps{pati}(:,:,:,cvali);
            if any(size(prediction.data.outmaps{pati})~=size(rsi_predicted(:,:,:,:,pati)))
                rsi_predicted(:,:,:,cvali,pati) = imresize3(f(cvali)*prediction.data.outmaps{pati}(:,:,:,cvali), size(rsi_predicted(:,:,:,cvali,pati)));

            else
                rsi_predicted(:,:,:,cvali,pati) = f(cvali)*prediction.data.outmaps{pati}(:,:,:,cvali);
            end
            calcDWI(:,:,:,bvali,pati) = calcDWI(:,:,:,bvali,pati) + rsi_predicted(:,:,:,cvali,pati).*exp(-bvals(bvali).*params.ModelADCs(cvali));
        end
    end

    idx = find(bvals == 0);
    tmp = calcDWI(:,:,:,idx,pati);
    mb0(pati)= median(tmp(ProstateSegmentation(:,:,:,pati)>0.5),'all');
end

%     scalingFactor = regress(mb0_ref,mb0);
%     v1 = mb0_ref;
%     v2 = mb0;
%
%     defvec = isfinite(v1+v2);
%     b = pinv(v2(defvec))*v1(defvec);
%
%     v1_hat = b*v2;
%
%     disp([b mean((v1(defvec)-v2(defvec)).^2)/mean(v1(defvec).^2) mean((v1(defvec)-v1_hat(defvec)).^2)/mean(v1(defvec).^2)])
%
%     figure; plot(v1,v2,'*',[0 500],[0 500],'k:','LineWidth',2); axis equal; lim = [0 max([xlim ylim])]; xlim(lim); ylim(lim); h = [xlabel('mb0 TEmin') ylabel(sprintf('%s original','mb0 TE90')) title('mb0 calibration')]; set(h,'FontSize',14,'FontWeight','bold','Interpreter','none');
%
%     figure; plot(v1,v1_hat,'*',[0 500],[0 500],'k:','LineWidth',2); axis equal; lim = [0 max([xlim ylim])]; xlim(lim); ylim(lim); h = [xlabel('mb0 TEmin') ylabel(sprintf('%s corrected', 'mb0 TEmin*')) title('mb0 calibration')]; set(h,'FontSize',14,'FontWeight','bold','Interpreter','none');
%

end

