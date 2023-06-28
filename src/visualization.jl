export plotcamera

function plotcamera(Ci, l, parent; col=[0, 0, 1], plotCamPath=false)
    # from https://github.com/peterkovesi/ImageProjectiveGeometry.jl/blob/master/src/projective.jl#L2149
    scene = Scene(parent)
    if isa(Ci, Array)
        C = Ci
    else
        C = Vector(1)
        C[1] = Ci
    end

    for i = 1:length(C)

        if C[i].rows == 0 || C[i].cols == 0
            @warn("Camera rows and cols not specified")
            continue
        end

        f = C[i].fx  # Use fx as the focal length

        if i > 1 && plotCamPath
            lines!(scene, [C[i-1].P[1], C[i].P[1]],
                [C[i-1].P(2), C[i].P[2]],
                [C[i-1].P(3), C[i].P[3]])
        end

        # Construct transform from camera coordinates to world coords
        Tw_c = [C[i].Rc_w' C[i].P
            0 0 0 1]

        # Generate the 4 viewing rays that emanate from the principal point and
        # pass through the corners of the image.
        ray = zeros(3, 4)
        ray[:, 1] = [-C[i].cols / 2, -C[i].rows / 2, f]
        ray[:, 2] = [C[i].cols / 2, -C[i].rows / 2, f]
        ray[:, 3] = [C[i].cols / 2, C[i].rows / 2, f]
        ray[:, 4] = [-C[i].cols / 2, C[i].rows / 2, f]

        # Scale rays to distance l from the focal plane and make homogeneous
        ray = makehomogeneous(ray * l / f)
        ray = Tw_c * ray                 # Transform to world coords

        for n = 1:4             # Draw the rays
            lines!(scene, [C[i].P[1], ray[1, n]],
                [C[i].P[2], ray[2, n]],
                [C[i].P[3], ray[3, n]],
                color=col)
        end

        # Draw rectangle joining ends of rays
        lines!(scene, [ray[1, 1], ray[1, 2], ray[1, 3], ray[1, 4], ray[1, 1]],
            [ray[2, 1], ray[2, 2], ray[2, 3], ray[2, 4], ray[2, 1]],
            [ray[3, 1], ray[3, 2], ray[3, 3], ray[3, 4], ray[3, 1]],
            color=col)

        # Draw and label axes
        X = Tw_c[1:3, 1] * l .+ C[i].P
        Y = Tw_c[1:3, 2] * l .+ C[i].P
        Z = Tw_c[1:3, 3] * l .+ C[i].P

        lines!(scene, [C[i].P[1], X[1, 1]], [C[i].P[2], X[2, 1]], [C[i].P[3], X[3, 1]],
            color=col)
        lines!(scene, [C[i].P[1], Y[1, 1]], [C[i].P[2], Y[2, 1]], [C[i].P[3], Y[3, 1]],
            color=col)
        #    plot3D([C[i].P[1], Z(1,1)], [C[i].P[2], Z(2,1)], [C[i].P[3], Z(3,1)],...
        #           color=col)
        # text3D(X[1], X[2], X[3], "X", color=col)
        # text3D(Y[1], Y[2], Y[3], "Y", color=col)

    end
    scene
end
