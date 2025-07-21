function graph_cam(pose_cam,type)
    axes =  [0, 0, 0, 1; 100, 0, 0, 1; 0, 100, 0, 1; 0, 0, 100, 1]';
    p2d = pose_cam(1:3, 1:4) * axes;
    line([p2d(1,1) p2d(1,2)],[p2d(2,1) p2d(2,2)],[p2d(3,1) p2d(3,2)],'Color','red');
    hold on
    line([p2d(1,1) p2d(1,3)],[p2d(2,1) p2d(2,3)],[p2d(3,1) p2d(3,3)],'Color','green');
    hold on
    line([p2d(1,1) p2d(1,4)],[p2d(2,1) p2d(2,4)],[p2d(3,1) p2d(3,4)],'Color','blue');
    if type==1
        text(p2d(1,1),p2d(2,1),p2d(3,1),'VIS-Cam')
    elseif type==2
        text(p2d(1,1),p2d(2,1),p2d(3,1),'IR-Cam')
    elseif type==3
        text(p2d(1,1),p2d(2,1),p2d(3,1),'Proj')
    end
end