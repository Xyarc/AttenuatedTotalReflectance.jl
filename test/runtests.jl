using AttenuatedTotalReflectance
using LinearAlgebra
using StaticArrays
using Test

@testset "AttenuatedTotalReflectance.jl" begin
    
    @testset "Utility Functions" begin
        
        # Testing complex_n
        @test complex_n(1.5, 0.5) === 1.5 + 0.5im
        @test complex_n(1.0, 0.0) isa ComplexF64

        @testset "Snell's Law" begin

            @testset "Standard Refraction (Real Angles)" begin
                # Air to Glass: n_in=1.0, n_out=1.5, theta=30°
                # Use deg2rad because sin() expects radians
                theta_in = deg2rad(30.0)
                n_in = 1.0
                n_out = 1.5

                result = snells_law(theta_in, n_in, n_out)

                @test result isa Complex
                @test real(result) ≈ asin(1.0 / 1.5 * sin(theta_in))
                @test imag(result) ≈ 0.0 atol = 1e-15
            end

            @testset "Total Internal Reflection (TIR)" begin
                # Glass to Air: n_in=1.5, n_out=1.0, theta=60° (Critical angle is ~41.8°)
                theta_in = deg2rad(60.0)
                n_in = 1.5
                n_out = 1.0

                result = snells_law(theta_in, n_in, n_out)

                # In TIR, the real part of the angle should be π/2 (90°)
                @test real(result) ≈ π / 2
                # The imaginary part should be non-zero (representing the evanescent component)
                @test imag(result) > 0
            end

            @testset "Normal Incidence" begin
                # theta = 0 should always return 0
                @test snells_law(0.0, 1.5, 1.0) ≈ 0.0
            end

            @testset "Type Stability" begin
                # Testing that the compiler can infer the Complex return type
                @test (@inferred snells_law(0.5, 1.0, 1.5)) isa Complex
            end
            @testset "Critical Angle Transition" begin
                n_in = 1.5   # Denser (Glass)
                n_out = 1.0  # Rarer (Air)

                # 1. Calculate the theoretical critical angle
                theta_c = asin(n_out / n_in) # ~0.7297 rad or 41.81°

                @testset "Exactly at Critical Angle" begin
                    result = snells_law(theta_c, n_in, n_out)

                    # Real part must be 90 degrees
                    @test real(result) ≈ π / 2
                    # Imaginary part should be zero (or effectively zero)
                    @test imag(result) ≈ 0.0 atol = 1e-15
                end

                @testset "Just Past Critical Angle (TIR)" begin
                    # 0.001 radians beyond critical
                    theta_tip = theta_c + 0.001
                    result = snells_law(theta_tip, n_in, n_out)

                    # Real part stays locked at 90 degrees (π/2)
                    @test real(result) ≈ π / 2
                    # Imaginary part must now be non-zero
                    @test imag(result) > 0.0
                end

                @testset "Just Before Critical Angle" begin
                    # 0.001 radians before critical
                    theta_safe = theta_c - 0.001
                    result = snells_law(theta_safe, n_in, n_out)

                    # Should be purely real and slightly less than π/2
                    @test real(result) < π / 2
                    @test imag(result) ≈ 0.0 atol = 1e-15
                end
            end
            @testset "Edge Cases" begin
                # Equal indices (theta_in should equal theta_out)
                @test real(snells_law(0.8, 1.33, 1.33)) ≈ 0.8

                # Test with Float32 to ensure it doesn't force a promotion to Float64 if not needed
                # (Though sin/asin usually promote to Float64 in Julia)
                @test snells_law(0.5f0, 1.0f0, 1.5f0) isa Complex
            end
        end

        @testset "Permittivity to NK Conversion" begin

            @testset "Mathematical Correctness" begin
                # Case: Vacuum/Air (e1=1, e2=0) -> n=1, k=0
                n, k = epsilon_to_nk(1.0, 0.0)
                @test n ≈ 1.0
                @test k ≈ 0.0

                # Case: Perfect Metal/Plasma at resonance (e1=0, e2=1)
                # n = sqrt(0.5 * (1 + 0)) = sqrt(0.5)
                # k = sqrt(0.5 * (1 - 0)) = sqrt(0.5)
                n, k = epsilon_to_nk(0.0, 1.0)
                @test n ≈ sqrt(0.5)
                @test k ≈ sqrt(0.5)

                # Test consistency: (n + ik)^2 should equal (e1 + ie2)
                e1, e2 = -10.5, 2.3
                n, k = epsilon_to_nk(e1, e2)
                @test n^2 - k^2 ≈ e1
                @test 2 * n * k ≈ e2
            end

            @testset "Type Stability" begin
                # Ensure it always returns Float64 even with weird inputs
                @inferred epsilon_to_nk(10.0, 5.0)
                n, k = epsilon_to_nk(1.0, 1.0)
                @test n isa Float64
                @test k isa Float64
            end

            @testset "Edge Cases" begin
                # Large values (High-index materials)
                @test all(isfinite.(epsilon_to_nk(1e6, 1e6)))

                # Zero permittivity
                @test epsilon_to_nk(0.0, 0.0) == (0.0, 0.0)
            end
        end

        @testset "NK to Permittivity Conversion" begin
    
            @testset "Mathematical Correctness" begin
            # Case: Vacuum (n=1, k=0) -> e1=1, e2=0
            e1, e2 = nk_to_epsilon(1.0, 0.0)
            @test e1 == 1.0
            @test e2 == 0.0

            # Case: Simple Dielectric (n=1.5, k=0) -> e1=2.25, e2=0
            e1, e2 = nk_to_epsilon(1.5, 0.0)
            @test e1 ≈ 2.25
            @test e2 == 0.0

            # Case: Absorbing Material (n=1, k=1) -> e1=0, e2=2
            e1, e2 = nk_to_epsilon(1.0, 1.0)
            @test e1 == 0.0
            @test e2 == 2.0

                @testset "Type Stability & Precision" begin

                    @inferred nk_to_epsilon(1.44, 0.001)

                    # Check that it handles zero correctly
                    @test nk_to_epsilon(0.0, 0.0) == (0.0, 0.0)

                    # Large values test
                    e1, e2 = nk_to_epsilon(1e4, 1e4)
                    @test isfinite(e1) && isfinite(e2)
                end
            end
        end


        @testset "Optical Constant Round-trip Tests" begin
            # Test that nk_to_epsilon(epsilon_to_nk(e1, e2)) returns (e1, e2)
            original_e1, original_e2 = -10.5, 2.3
            n, k = epsilon_to_nk(original_e1, original_e2)
            final_e1, final_e2 = nk_to_epsilon(n, k)

            @test final_e1 ≈ original_e1
            @test final_e2 ≈ original_e2
        end
    

    end

    @testset "Layer Construction" begin
        
        # Testing Layer Struct
        l = layer("Silver", 0.051585, 3.9046, 20e-9)

        @test l.material == "Silver"
        @test l.thickness == 20e-9
        @test l.n ≈ 0.051585
        @test l.k ≈ 3.9046

        # Test for error on negative thickness
        @test_throws ArgumentError layer("Bad Layer", 0.0, 1.0, -10.0)


        # Define some dummy materials
        mats = [
            (material="Glass", n=1.5, k=0.0, thickness=5e-6),
            (material="Gold", n=0.18, k=3.0, thickness=50e-9),
            (material="Air", n=1.0, k=0.0, thickness=5e-6),
        ]

            stack = material_stack(mats)

        @testset "Structure Integrity" begin
            # Substrate + 1 metal + Last Dielectric Layer + PML = 4
            @test length(stack) == 4
            @test stack[1].material == "Substrate: Glass"
            @test stack[end].material == "PML"
            
        end

        @testset "Data Pass-through" begin
            # Verify the second material (Gold) is in the correct slot
            @test stack[end-1].material == "Air"
            @test stack[2].n == 0.18
            @test stack[2].thickness ≈ 50e-9
        end

    end


    @testset "Core Physics & TMM" begin
        
        @testset "Target Layer Identification" begin
            # Setup a standard 4-layer stack
            stack = [
                layer("Prism", 1.5, 0.0, 1e-6),
                layer("Gold", 0.18, 3.0, 50e-9,),
                layer("FlowCell", 1.33, 0.0, 5e-6),
                layer("PML", 1.33 + 1e-5, 0.0, 1e-6)
            ]

            @testset "Manual Selection" begin
                # Specifically target layer above the metal layer (Index 2)
                idx, n = target_layer(stack, 0.0, 2)
                @test idx == 3
                @test n ≈ 1.33 + 0im
            end

            @testset "No Layer Specified" begin
                # If 0 is passed, it should pick length(stack) - 1 (The FlowCell/Index 3)
                
                @test_throws ArgumentError target_layer(stack, 0.0, 0)
            end

            @testset "Inference" begin
                # Crucial for the speed of your ATR loops
                @test (@inferred target_layer(stack, 0.0, 2)) isa Tuple{Int,ComplexF64}
            end
        end

        # Mock Data for testing
        n_air = 1.0 + 0.0im
        n_glass = 1.5 + 0.0im
        wavelength = 633.0
        theta_30 = deg2rad(30.0)
        theta_glass = snells_law(theta_30, n_air, n_glass)

        @testset "Fresnel Coefficients" begin
            # Test normal incidence (S and P should be identical at 0 degrees)
            rs = refl_coeff_S(0.0, n_air, 0.0, n_glass)
            rp = refl_coeff_P(0.0, n_air, 0.0, n_glass)
            @test rs ≈ rp 
            @test abs2(rs) ≈ ((1.5 - 1.0) / (1.5 + 1.0))^2

        end

        @testset "Transfer Matrix Properties" begin
            # Simple 3-layer stack: Air / 50nm Gold / Air
            # Note: Using your Layer struct
            stack = [
                layer("Air", 0.0, 1.0, 0.0),
                layer("Gold", 50.0, 0.18, 3.0),
                layer("Air", 0.0, 1.0, 0.0)
            ]

            C, D0, Dj, t_tot, t0, tj = compute_transfer_coefficients(stack, theta_30, wavelength)

            @test C isa SMatrix{2,2,ComplexF64}
            # Determinant of a T-matrix for non-magnetic media should be the ratio of admittance
            # For simple cases, we just check it's not zero and is type-stable
            @test !iszero(det(C))
        end

        @testset "ATR Simulation (Angular)" begin
            # Kretschmann Configuration: Prism (1.5) / Gold (50nm) / Air (1.0)
            # Resonance for Gold at 633nm is usually around 43-45 degrees
            prism_stack = [
                layer("Prism", 1.51, 0.0,  1e-5),
                layer("Silver",  0.056206, 4.2776, 47e-9), # typical nk for Silver at 633 nm
                layer("Air", 1.0, 0.0, 5e-6),
                layer("PML", 1.0+1e-5, 0.0, 1e-6)
            ]

            angles = deg2rad.(collect(35.0:0.01:55.0))
            R, T, fx, fy, fz, fp = angular_ATR(prism_stack, angles, 633.0e-9)

            @test length(R) == length(angles)
            @test all(0 .<= R .<= 1.0)
            @test all(0 .<= T .<= 1.0)

            # Check that field enhancement fp is non-negative
            @test all(fp .>= 0.0)

            # Physical check: At 30 degrees (below critical angle ~41.8), T should be > 0
            @test T[1] > 0.001
            # Physical check: At 50 degrees (TIR), T should be near 0
            @test T[end] < 1e-5
        end

        @testset "Inference & Performance" begin
            stack = [layer("Test", 10.0, 1.5, 0.0)]
            # Ensure no type-instability in the core loop
            @test (@inferred compute_transfer_coefficients(stack, 0.5, 633.0)) isa Tuple
        end
    end
end
