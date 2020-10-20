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

(define (show-2pt traj)
  (define ld00 (load-two-point traj 0 0))
  (define ld01 (load-two-point traj 0 1))
  (define ld10 (load-two-point traj 1 0))
  (define ld11 (load-two-point traj 1 1))
  (print-datatable
    (vector-map
      (lambda (n)
        (vector
          n
          (real-part
            (vector-head (lat-data-ref ld00 (vector n 15 15))))
          (real-part
            (vector-head (lat-data-ref ld01 (vector n 15 15))))
          (real-part
            (vector-head (lat-data-ref ld10 (vector n 15 15))))
          (real-part
            (vector-head (lat-data-ref ld11 (vector n 15 15))))
          ))
      (list->vector (iota (t-size))))))

(define (show traj)
  (show-3pt-type traj 1 1 0)
  ; (show-2pt traj)
  )

(show 1900)
; (show 2260)
; (show 2270)

