#!/usr/bin/env scheme-script

(import (cslib all))

(define (get-a-iwasaki beta)
  (let* ([c1 -2.1281]
         [c2 -1.0056]
         [c3 0.6041]
         [b-3 (- beta 3)]
         [b-3^2 (sqr b-3)]
         [lna/r0 (+ c1 (* c2 b-3) (* c3 b-3^2))])
    (* 0.5 (exp lna/r0))))

(define (get-a-dbw2 beta)
  (let* ([d1 -1.6007]
         [d2 -2.3179]
         [d3 -0.8020]
         [d4 -19.8509]
         [b-1 (- beta 1)]
         [b-1^2 (sqr b-1)]
         [b-1^3 (* b-1 b-1^2)]
         [lna/r0 (+ d1 (* d2 b-1) (* d3 b-1^2) (* d4 b-1^3))])
    (* 0.5 (exp lna/r0))))

(print "hello world")

(print (get-a-iwasaki 2.31))

(print (get-a-dbw2 0.7796))

; (print (get-a-iwasaki 2.7124))
