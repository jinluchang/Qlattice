#!/usr/bin/env scheme-script

(import (cslib all))

; table format
; mass lx ly lz ts val-short val

(define table
  (list-sort
    (on < cadr)
    (map
      (lambda (s)
        (map string->number (cdr (words s))))
      (filter
        (lambda (s)
          (string-prefix? "RESULT:" s))
        (load-lines "log.txt")))))

(define table-0.14-t=l/2
  (filter
    (lambda (l)
      (and (= 0.14 (list-ref l 0))
           (= (* 2 (list-ref l 4)) (list-ref l 1))))
    table))

(define (get-val-m-l mass lsize)
  (list-nref
    (filter
      (lambda (l)
        (and
          (= mass (list-ref l 0))
          (= lsize (list-ref l 1))
          (= lsize (* 2 (list-ref l 4)))))
      table)
    0 6))

(define (filter-mass mass tab)
  (filter
    (lambda (l)
      (= mass (list-ref l 0)))
    tab))

(define (filter-l lsize tab)
  (filter
    (lambda (l)
      (= lsize (list-ref l 1)))
    tab))

(define (filter-ts ts tab)
  (filter
    (lambda (l)
      (= ts (list-ref l 4)))
    tab))

(define (get-val-short tab)
  (list-nref tab 0 5))

(define (get-val tab)
  (list-nref tab 0 6))

(define inf-val-0.14
  (get-val
    (filter-ts
      48
      (filter-mass
        0.14
        (filter-l 96 table)))))

(define (mk-dm-table tab inf-val)
  (define (f l)
    (pmatch l
      [(,mass ,lx ,ly ,lz ,ts ,val-short ,val)
       (list (* (/ mass 0.14) lx /fm/gev)
             (* 140 (/ val mass))
             (* 140 (/ (- val inf-val) mass)))]))
  (map f tab))

; (print table)
; (print table-0.14-t=l/2)

(plot-save
  "plots/size-mass.pdf.eps.png"
  "set logscale y"
  "set key tm"
  "set size 0.6,0.7"
  "set xlabel '$L$ (fm)'"
  ; "set ylabel 'finite volume error for E\\&M correction to $m_\\pi$ (MeV)'"
  "set ylabel 'FV correction to $\\Delta m_\\pi$ (MeV)'"
  (cons "size-dmass-tl2.txt" (table->datatable (mk-dm-table table-0.14-t=l/2 inf-val-0.14)))
  (mk-plot-line
    "plot [1.0:10.0] [:]"
    "'size-dmass-tl2.txt' u 1:3 w l t 'FV error with $t_s = L/2$'"
    )
  )

(plot-save
  "plots/ts-val.pdf.eps.png"
  "set key tm"
  "set size 0.6,0.9"
  "set xlabel '$t_s$ (fm)'"
  "set ylabel 'FV correction to $\\Delta m_\\pi$ (MeV)'"
  (cons "ts-val-32.txt"
        (table->datatable
          (map
            (lambda (l)
              (pmatch l
                [(,mass ,lx ,ly ,lz ,ts ,val-short ,val)
                 (list (* (/ mass 0.14) ts /fm/gev)
                       (* 140 (/ (- val-short inf-val-0.14) mass))
                       (* 140 (/ (- val inf-val-0.14) mass)))]))
          (filter
            (lambda (l)
              (and (= 0.14 (list-ref l 0))
                   (= 32 (list-ref l 1))))
            table))))
  (cons "ts-val-24.txt"
        (table->datatable
          (map
            (lambda (l)
              (pmatch l
                [(,mass ,lx ,ly ,lz ,ts ,val-short ,val)
                 (list (* (/ mass 0.14) ts /fm/gev)
                       (* 140 (/ (- val-short inf-val-0.14) mass))
                       (* 140 (/ (- val inf-val-0.14) mass)))]))
          (filter
            (lambda (l)
              (and (= 0.14 (list-ref l 0))
                   (= 24 (list-ref l 1))))
            table))))
  (mk-plot-line
    "plot [1.0:5.0] [-0.3:0.1]"
    "0 lc 1 not"
    "'ts-val-24.txt' u 1:3 w l lc 3 dt 1 t '$\\mathcal{I}^{\\phantom{(s)}}$ $L=4.7$ fm'"
    "'ts-val-24.txt' u 1:2 w l lc 3 dt 2 t '$\\mathcal{I}^{(s)}$ $L=4.7$ fm'"
    "'ts-val-32.txt' u 1:3 w l lc 2 dt 1 t '$\\mathcal{I}^{\\phantom{(s)}}$ $L=6.3$ fm'"
    "'ts-val-32.txt' u 1:2 w l lc 2 dt 2 t '$\\mathcal{I}^{(s)}$ $L=6.3$ fm'"
    )
  )

