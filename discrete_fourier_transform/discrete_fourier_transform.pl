use strict;
use warnings;
use Math::Complex;

sub naive_dft {
    my @x = @_;
    my $N = scalar(@x);
    my @X;

    for my $k (0 .. $N - 1) {
        my $sum_re = 0;
        my $sum_im = 0;
        for my $n (0 .. $N - 1) {
            my $angle = -2 * pi * $k * $n / $N;
            $sum_re += $x[$n] * cos($angle);
            $sum_im += $x[$n] * sin($angle);
        }
        push @X, Math::Complex->make($sum_re, $sum_im);
    }

    return @X;
}

# Example usage
my @x = (1, 2, 3, 4);
my @dft_result = naive_dft(@x);
print "DFT: ", join(", ", map { $_->Re . "+" . $_->Im . "i" } @dft_result), "\n";
