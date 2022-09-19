#!/usr/bin/env scheme-script

(import (cslib all))

(print "hello")

(define job-tag
  (make-parameter #f))

(define t-size
  (make-parameter #f))

(job-tag "24D")
(t-size 64)

(define (get-two-point traj)
  (format "analysis/lat-two-point/~a/results=~a" (job-tag) traj))

(define (get-two-point-name type1 type2)
  (format "two-point-wall-snk-sparse-corrected-~a-~a.lat" type1 type2))

(define (load-two-point traj type1 type2)
  (load-lat-data
    (filepath-append (get-two-point traj)
                     (get-two-point-name type1 type2))))

(define (get-meson-snk-src traj)
  (format "analysis/lat-meson-snk-src/~a/results=~a" (job-tag) traj))

(define (get-meson-snk-src-name type1 type2)
  (format "meson-snk-src-~a-~a.lat" type1 type2))

(define (load-meson-snk-src traj type1 type2)
  (load-lat-data
    (filepath-append (get-meson-snk-src traj)
                     (get-meson-snk-src-name type1 type2))))

(define (avg-meson-snk-src ld)
  (vector-map
    (lambda (tsep)
      (apply average
        (map (lambda (tsrc)
               (let ([tsnk (mod (+ tsep tsrc) (t-size))])
                 (vector-head (lat-data-ref ld (vector tsnk tsrc)))))
             (iota (t-size)))))
    (list->vector (iota (t-size)))))

(define (get-three-point traj)
  (format "analysis/lat-three-point/~a/results=~a" (job-tag) traj))

(define (get-three-point-name type1 type2 type3)
  (format "three-point-~a-~a-~a.lat" type1 type2 type3))

(define (load-three-point traj type1 type2 type3)
  (load-lat-data
    (filepath-append (get-three-point traj)
                     (get-three-point-name type1 type2 type3))))

(define (show-3pt-type traj type1 type2 type3)
  (define ld-2pt-1 (load-two-point traj type1 type3))
  (define ld-2pt-2 (load-two-point traj type2 type3))
  (define ld-3pt (load-three-point traj type1 type2 type3))
  (print-datatable
    (vector-map
      (lambda (tsep)
        (vector-append
          (vector tsep)
          (vector-map
            (lambda (top)
              (/ (sqrt (* (real-part (vector-head (lat-data-ref ld-2pt-1 (vector (* 2 top) 15 15))))
                          (real-part (vector-head (lat-data-ref ld-2pt-2 (vector (* 2 (- tsep top)) 15 15))))))
                 (real-part (vector-head (lat-data-ref ld-3pt (vector tsep top 8))))))
            (list->vector (make-seq 1 (- tsep 1))))))
      (list->vector (make-seq 2 (/ (t-size) 3))))))

(define (show-2pt-type traj type1 type2)
  (define ld (load-two-point traj type1 type2))
  (define v-mss (avg-meson-snk-src (load-meson-snk-src traj type1 type2)))
  (print-datatable
    (vector-map
      (lambda (n)
        (vector
          n
          (real-part
            (vector-head (lat-data-ref ld (vector n 15 15))))
          (real-part
            (vector-ref v-mss n))
          (/ (real-part
               (vector-ref v-mss n))
             (real-part
               (vector-head (lat-data-ref ld (vector n 15 15)))))
          ))
      (list->vector (iota (t-size))))))

(define (show traj)
  (show-2pt-type traj 0 0)
  ; (show-3pt-type traj 1 1 0)
  )

(show 1010)
; (show 1900)
; (show 2260)
; (show 2270)

