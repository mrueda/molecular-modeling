use strict;
use warnings;
use Math::Complex;

sub fft_recursive {
    my @x = @_;
    my $N = scalar(@x);

    # Base case: If only one element, return it
    return @x if $N == 1;

    # Input length must be a power of 2
    die "Input length must be a power of 2" if $N % 2 != 0;

    # Divide into even and odd parts
    my @even = fft_recursive(map { $x[$_] } grep { $_ % 2 == 0 } (0 .. $N - 1));
    my @odd  = fft_recursive(map { $x[$_] } grep { $_ % 2 == 1 } (0 .. $N - 1));

    # Combine using butterfly operations
    my @X;
    for my $k (0 .. $N / 2 - 1) {
        my $t = $odd[$k] * Math::Complex->make(cos(-2 * pi * $k / $N), sin(-2 * pi * $k / $N));
        $X[$k]          = $even[$k] + $t;
        $X[$k + $N / 2] = $even[$k] - $t;
    }

    return @X;
}

# Example usage
my @x = (1, 2, 3, 4);
#my @x = (1, 2, 3, 4, 0, 0, 0, 0);  # Zero-padding to make length a power of 2
my @fft_result = fft_recursive(@x);
print "FFT: ", join(", ", map { $_->Re . "+" . $_->Im . "i" } @fft_result), "\n";
