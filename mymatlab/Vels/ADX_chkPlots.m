hSum=figure;
pcolor(ADX.yday_DN,ADX.depth,(ADX.ec1+ADX.ec2+ADX.ec3+ADX.ec4)/4),axis ij,shading flat
hold on
plot(ADX.yday_BT,ADX.bottomBT,'k-',ADX.yday_BT,ADX.bottomBT,'k.')
hW=figure; clf
hBm=figure; clf
hVe=figure; clf
%%% Set bad/error parameters
EcAR = [175,25]; % echo>A,dist>R = bad
difPtRl_Thresh = 4; % badflag pings where pitch/roll changed by more, ping-to-ping

ix = find( abs(diff(ADX.pitch_DN)) > difPtRl_Thresh | ...
    abs(diff(ADX.roll_DN)) > difPtRl_Thresh | ...
    abs(diff(ADX.pitch_UP)) > difPtRl_Thresh | ...
    abs(diff(ADX.roll_UP)) > difPtRl_Thresh );
figure(hSum)
if ~isempty(ix) % pings ignored because of pitching
    plot(ADX.yday_DN(ix+1),ADX.pr_DN(ix+1)*100, 'kx')
end
i = 1; GetOut = 0;
while i<length(ADX.yday_DN)+1
    ig = find(~isnan(ADX.w(:,i)));
    if length(ig)>2
        figure(hSum)
        plot(ADX.yday_DN(i),ADX.pr_DN(i)*100,'ro')
        rBot = ADX.bottomBT(i)-ADX.pr_BT(i)*100;
        ieR = find(ADX.status(:,i) > 0 & ADX.status(:,i) < 1024); % some status bits set
        ieE = find(ADX.status(:,i) == 32); % echo amp is only one set
        ieC = find(ADX.status(:,i) == 1024); % rejected as too few contiguous bins
%         ieE = find(~isnan(ADX.bm1(:,i)) & ...
%             max([ADX.bm1(:,i) ADX.bm2(:,i) ADX.bm3(:,i) ADX.bm4(:,i)], [], 2) > EcAR(1) & ...
%             (ADX.depth < ADX.pr_UP(i)*100-EcAR(2) | ADX.depth > ADX.pr_DN(i)*100+EcAR(2)) );
        figure(hW),clf
        zb = ADX.depth(ig(1));ze=ADX.depth(ig(end));
        plot(ADX.w(ig(1):ig(end),i),ADX.depth(ig(1):ig(end)),'b-',...
            ADX.w(ig(1):ig(end),i),ADX.depth(ig(1):ig(end)),'b.')
        hold on, axis ij
        plot(ADX.w(ieR,i),ADX.depth(ieR),'r.',ADX.w(ieE,i),ADX.depth(ieE),'mo', ...
            ADX.w(ieC,i),ADX.depth(ieC),'r^')
        %plot(ADX.werr(ig(1):ig(end),i),ADX.depth(ig(1):ig(end)),'y-'), axis ij
        plot(ADX.wBT(i),ADX.pr_BT(i)*100,'ko')
        plot([ADX.w_UP(i) ADX.w_UP(i)], [zb ze], 'g-',[ADX.w_DN(i) ADX.w_DN(i)], [zb ze], 'r-')
        iws = [max(1,i-2):min(length(ADX.yday_DN),i+2)];
        plot(ADX.w_UP(iws),ADX.pr_UP(iws)*100,'g-',ADX.w_DN(iws),ADX.pr_DN(iws)*100,'r-')
        title(['col.' num2str(i) ' ydDN=' num2str(mod(ADX.yday_DN(i),1))]), grid on
        figure(hBm),clf
        ig = find(~isnan(ADX.ec1(:,i)));
        plot(ADX.ec1(ig(1):ig(end),i),ADX.depth(ig(1):ig(end)),'k-'), hold on,axis ij
        plot(ADX.ec1(ieR,i),ADX.depth(ieR),'k.',ADX.ec1(ieE,i),ADX.depth(ieE),'ko')
        plot(ADX.ec2(ig(1):ig(end),i),ADX.depth(ig(1):ig(end)),'r-')
        plot(ADX.ec2(ieR,i),ADX.depth(ieR),'r.',ADX.ec2(ieE,i),ADX.depth(ieE),'ro')
        plot(ADX.ec3(ig(1):ig(end),i),ADX.depth(ig(1):ig(end)),'b-')
        plot(ADX.ec3(ieR,i),ADX.depth(ieR),'b.',ADX.ec3(ieE,i),ADX.depth(ieE),'bo')
        plot(ADX.ec4(ig(1):ig(end),i),ADX.depth(ig(1):ig(end)),'g-')
        plot(ADX.ec4(ieR,i),ADX.depth(ieR),'g.',ADX.ec4(ieE,i),ADX.depth(ieE),'go')
        grid on
        figure(hVe),clf
        ig = find(~isnan(ADX.u(:,i)+ADX.v(:,i))); % plot -v, then (u,-v)>0, ord(.71*fwd vel)
        plot(ADX.u(ig(1):ig(end),i),ADX.depth(ig(1):ig(end)),'g-'), hold on,axis ij
        plot(-ADX.v(ig(1):ig(end),i),ADX.depth(ig(1):ig(end)),'r-')
        plot([ADX.uBT(i) ADX.uBT(i)], [zb ze], 'c:',[-ADX.vBT(i) -ADX.vBT(i)], [zb ze], 'm:')
        plot(ADX.u(ieR,i),ADX.depth(ieR),'g.',ADX.u(ieE,i),ADX.depth(ieE),'co')
        plot(-ADX.v(ieR,i),ADX.depth(ieR),'r.',-ADX.v(ieE,i),ADX.depth(ieE),'mo')
        grid on

        if ~mod(i,40)
            x=input('cont (>0 to set GetOut=1): ');
            if x>0
                GetOut = 0
                keyboard
            end
        else
            pause
        end
    end
    i = i+1;
    if GetOut
        break
    end
end