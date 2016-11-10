yd_b = 87.523; yd_e = 400;
e_b = max(ADDN.ens_no(1), ADUP.ens_no(1)); e_b=560;
e_e = min(ADDN.ens_no(end), ADUP.ens_no(end));

id = find(ADDN.yday>=yd_b & ADDN.yday<=yd_e & ...
    ADDN.ens_no>=e_b & ADDN.ens_no<=e_e);
ADX.ens_no = ADDN.ens_no(id);
ADX.yday = ADDN.yday(id);
iu = find(ADUP.ens_no>=ADX.ens_no(1)& ADUP.ens_no<=ADX.ens_no(end));

tvars = {'bottomBT','VbtE','VbtN','VbtZ','SWIMS_pitch','SWIMS_roll'};
for i=1:length(tvars)
    eval(['ADX.' tvars{i} '=ADDN.' tvars{i} '(id);']);
end

avars = {'v1_bm','v2_bm','v3_bm','v4_bm', ...
        'ec1_bm','ec2_bm','ec3_bm','ec4_bm', ...
        'cor1_bm','cor2_bm','cor3_bm','cor4_bm'};
icX = [50:-1:1];
ADX.z_off = [ADUP.z_adcp(icX) 0 ADDN.z_adcp];
zrow = NaN*iu;
for i=1:length(avars)
    eval(['ADX.' avars{i} '=[ADUP.' avars{i} '(icX,iu); zrow; ADDN.' avars{i} '(:,id)];']);
end

return

figure,pcolor(ADX.yday,ADX.z_off,ADX.cor2_bm),shading flat,caxis([0 140]),colormap(jet),colorbar,axis ij
hold on
plot(SWraw.SWIMStime,-100*SWraw.Pr,'k-')
hold on
hold on
plot(SWraw.SWIMStime,-100*SWraw.Pr,'k-')
plot(ADX.yday,ADX.bottomBT,'k-')
set(gca,'ylim',[-150 220])
plot(SWraw.SWIMStime,-100*SWraw.Pr+120,'w-')
plot(SWraw.SWIMStime,-100*SWraw.Pr+150,'g-')
figure,pcolor(ADX.yday,ADX.z_off,ADX.ec2_bm),shading flat,colormap(jet),colorbar,axis ij
plot(SWraw.SWIMStime,-100*SWraw.Pr+150,'g-')
clf
close
figure,pcolor(ADX.yday,ADX.z_off,ADX.ec2_bm),shading flat,colormap(jet),colorbar,axis ij
hold on
plot(SWraw.SWIMStime,-100*SWraw.Pr+150,'g-')
plot(ADX.yday,ADX.bottomBT,'k-')
plot(SWraw.SWIMStime,-100*SWraw.Pr,'k-')
set(gca,'ylim',[-150 220])