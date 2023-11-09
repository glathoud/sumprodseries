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

enum ATKIN_UPDATE = false;

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

  static if (!ATKIN_UPDATE)
    scope big_primes = get_big_primes();
  
  foreach (i; 0..n_max)
    {
      auto x = r.front;

      static if (!ATKIN_UPDATE)
        if (x > big_primes[ $-1 ])
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

  atkin_update_big_primes( big_primes, max( x, big_primes[ $-1 ] + 1000 ) );
                           
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



void atkin_update_big_primes( ref BigInt[] big_primes, BigInt limitSqrt )
{
  /* based on:
     https://en.wikipedia.org/wiki/Sieve_of_Atkin
     https://gist.github.com/kramtark/277c0657530572e92c48
     (see copy below)
   */

  // xxx todo: consider taking `limit` as input and computing `limitSqrt`. Less overkill.
  
  scope auto limit = limitSqrt * limitSqrt;

  // Move from `big_primes` representation to `sieve` representation
  scope bool[BigInt] sieve;
  foreach (bp; big_primes)  sieve[ bp ] = true;
  
  scope BigInt n, x, xx, y, yy, i;

  scope immutable bp_init_max = big_primes[$-1];
  
  for (x = 2 + bp_init_max; x <= limitSqrt; ++x)
    {
      xx = x*x;
      for (y = 1; y <= limitSqrt; ++y) {
        yy = y*y;
        if (xx + yy >= limit) {
          break;
         }
         // first quadratic using m = 12 and r in R1 = {r : 1, 5}
         n = (4 * xx) + (yy);
         if (n <= limit && (n % 12 == 1 || n % 12 == 5)) {
           sieve[n] = !sieve[n];
         }
         // second quadratic using m = 12 and r in R2 = {r : 7}
         n = (3 * xx) + (yy);
         if (n <= limit && (n % 12 == 7)) {
           sieve[n] = !sieve[n];
         }
         // third quadratic using m = 12 and r in R3 = {r : 11}
         n = (3 * xx) - (yy);
         if (x > y && n <= limit && (n % 12 == 11)) {
           sieve[n] = !sieve[n];
         }
       }
    }

  // false each primes multiples
  for (n = 5; n <= limitSqrt; n++) {
    if (sieve[n]) {
      x = n * n;
      for (i = x; i <= limit; i += x) {
        sieve[i] = false;
      }
    }
  }
  
  // Actual update: move back from `sieve` representation to `big_primes` representation
  scope auto bp_app = appender!(BigInt[]);
  
  bp_app.put( big_primes );
  
  foreach (ref bp; sieve.keys.dup.sort)
    {
      if (bp > bp_init_max  &&  sieve[ bp ])
        bp_app.put( bp );
    }
  
  big_primes = bp_app.data;
}

/* original JS implementation

 function sieveOfAtkin(limit){
   var limitSqrt = Math.sqrt(limit);
   var sieve = [];
   var n;

   //prime start from 2, and 3
   sieve[2] = true;
   sieve[3] = true;

   for (var x = 1; x <= limitSqrt; x++) {
       var xx = x*x;
       for (var y = 1; y <= limitSqrt; y++) {
           var yy = y*y;
           if (xx + yy >= limit) {
             break;
           }
           // first quadratic using m = 12 and r in R1 = {r : 1, 5}
           n = (4 * xx) + (yy);
           if (n <= limit && (n % 12 == 1 || n % 12 == 5)) {
               sieve[n] = !sieve[n];
           }
           // second quadratic using m = 12 and r in R2 = {r : 7}
           n = (3 * xx) + (yy);
           if (n <= limit && (n % 12 == 7)) {
               sieve[n] = !sieve[n];
           }
           // third quadratic using m = 12 and r in R3 = {r : 11}
           n = (3 * xx) - (yy);
           if (x > y && n <= limit && (n % 12 == 11)) {
               sieve[n] = !sieve[n];
           }
       }
   }

   // false each primes multiples
   for (n = 5; n <= limitSqrt; n++) {
       if (sieve[n]) {
           x = n * n;
           for (var i = x; i <= limit; i += x) {
               sieve[i] = false;
           }
       }
   }

   //primes values are the one which sieve[x] = true
   return sieve;
}

primes = sieveOfAtkin(5000);
*/

