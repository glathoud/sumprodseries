#!../../software/ldc2_impl/bin/rdmd --shebang  -debug -g -gs -gf -link-defaultlib-debug
     // #!../ldc2_impl/bin/rdmd --shebang  -g -gs -gf -link-defaultlib-debug
     // #!../ldc2_impl/bin/rdmd --shebang  -debug -g -gs -gf -link-defaultlib-debug
     // (-debug compiles faster, plus lets check a few more things)

import primes_4m;
import d_glat.lib_search_bisection;
import std.algorithm;
import std.array;
import std.bigint;
import std.conv;
import std.format;
import std.math;
import std.range;
import std.stdio;

void main()
{
  writeln("xxx"); stdout.flush;
  
  display_some_about( Fibo() );
  
  display_some_about( SPS!"abc"() );
  display_some_about( SPS!"cba"() );
  display_some_about( SPS!"bac"() );
}

struct Fibo
{
  BigInt a = 0, b = 1;

  string getName() { return "Fibo"; }
  bool empty() { return false; }

  auto front() { return a; }
  void popFront() { auto new_b = a+b; a = b; b = new_b; }
};

struct SPS(string spec)
{
  string getName() { return format( "SPS!%s", spec ); }

  BigInt a = 1, b = 1, c = 1;
  
  bool empty() { return false; }
  auto front() { return a; }

  private bool _checked;
  void popFront()
  {
    if (!_checked)
      {
        assert( spec.length == 3 );
        scope tmp = spec.dup.to!(dchar[]);
        tmp.sort;
        assert( tmp == "abc" );
        _checked = true;
      }
    auto new_c = mixin(format("%s+%s*%s",spec[ 0 ], spec[ 1 ], spec[ 2 ]));
    a = b; b = c; c = new_c;
  }
};


void display_some_about(R)( R r, size_t n_max = 1000, size_t n_digit_max = 60 )
{
  display_some_about!R( r.getName(), r, n_max, n_digit_max );
}

void display_some_about(R)( in string name, R r, size_t n_max = 1000, size_t n_digit_max = 60 )
{
  writeln; writeln( name ); writeln; stdout.flush;
  
  immutable fmt_i = "%."~to!string( ceil( log(cast(double)( n_max )) / log(10.0) ) )~"d";
  immutable fmt_x = "%"~to!string(n_digit_max)~"s";

  auto xxx_big_primes = get_big_primes();
  
  foreach (i; 0..n_max)
    {
      auto x = r.front;

      if (x > xxx_big_primes[ $-1 ]) // xxx for now - because our implementation is too slow
        break;       

      immutable xs = format( "%d", x ); 
      
      if (xs.length > n_digit_max)
        break;

      immutable is_prime = get_is_prime( x );
      immutable  d_prev  = get_d_prev( x );
      immutable  d_next  = get_d_next( x );
      
      writefln( "%s#"~fmt_i~": "~fmt_x~"  is_prime: %d  d_prev:%+3d  d_next:%+3d", name, i, xs, is_prime ? 1 : 0, d_prev, d_next );
      
      r.popFront;
    }
}

bool get_is_prime( in BigInt x )
{
  auto big_primes = get_big_primes();

  if (x <= big_primes[ $-1 ])
    return -1 < search_bisection_exact( big_primes, x );

  /* search for more prime numbers

     probably not efficient at all, but results cached into the `big_primes` array

     xxx todo: replace/combine with Atkin sieve 
  */
  while (x > big_primes[ $-1 ])
    {
      scope auto y = big_primes[ $-1 ]+2;
      scope BigInt q, r;
      while (y <= x)
        {
          scope y_half = y / 2;
          bool  is_y_prime = true;
          foreach (bp; big_primes)
            {
              if (bp > y_half)
                break;

              divMod( y, bp, q, r );

              writeln( "xxx y ", y , "   bp ", bp,  "  q ", q, "  r ", r  );
              

              if (r == 0)
                {
                  is_y_prime = false;
                  break;
                }
            }

          if (is_y_prime)
            {
              writeln( "xxx y is prime!   ", y );
              big_primes ~= y;
            }
          else
            {
              writeln( "xxx y is NOT prime!    ", y );
            }

          y += 2;
        }
    }

  return get_is_prime( x );
}

BigInt get_d_prev( in BigInt x ) { return get_prev( x ) - x; }
BigInt get_d_next( in BigInt x ) { return get_next( x ) - x; }

BigInt get_prev( in BigInt x0 )
{
  // xxx could optimize a bit, when x0 already in _big_primes
  
  BigInt x = x0;
  while (0 < (--x))
    {
      if (get_is_prime( x ))
        return x;
    }

  return BigInt( 0 );
}

BigInt get_next( in BigInt x0 )
{
  // xxx could optimize a bit, when x0 already in _big_primes[0..$-1]

  BigInt x = x0;
  while (0 < (++x))
    {
      if (get_is_prime( x ))
        return x;
    }

  return BigInt( 0 );
}


private BigInt[] _big_primes;

BigInt[] get_big_primes()
{
  if (_big_primes.length == 0)
    _big_primes = primes.map!BigInt.array;

  return _big_primes;
}

