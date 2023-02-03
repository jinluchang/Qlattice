#!/usr/bin/env scheme-script

(import (cslib all))

(print "hello world")

(define min-traj
  (make-parameter #f))

(define max-traj-num
  (make-parameter #f))

(define data-selector
  (make-parameter #f))

(define data-tag
  (make-parameter #f))

(define data-avg
  (make-parameter #f))

(define block-size
  (make-parameter #f))

(define max-corr-length
  (make-parameter #f))

(define max-change-length
  (make-parameter #f))

; -----------------------------------------

(define (get-data-path)
  "results/gf_info")

(define (get-data-fns)
  (list-sort
    (on < car)
    (filter
      (lambda (p)
        (or (eq? #f (max-traj-num))
            (< (car p) (+ (min-traj) (max-traj-num)))))
      (map
        (lambda (p)
          (cons (string->number (car p)) (cdr p)))
        (paths->ppairs
          "traj="
          (glob-expand-at (get-data-path) "traj=*.lat")
          ".lat")))))

(define (subtract-avg data)
  (let ([avg (if (eq? #f (data-avg))
                 (apply average data)
                 (data-avg))])
    (map
      (lambda (x) (- x avg))
      data)))

(define (get-all-data)
  (map
    (lambda (p)
      ((data-selector)
       (load-lat-data (cdr p))))
    (get-data-fns)))

(define (get-data)
  (map
    (lambda (p)
      ((data-selector)
       (load-lat-data (cdr p))))
    (filter
      (lambda (p)
        (and
          (>= (car p) (min-traj))))
      (get-data-fns))))

(define (mk-corr xs)
  (let* ([n (length xs)]
         [np (- n (max-corr-length))]
         [xs1 (drop (max-corr-length) xs)])
    (lambda (i)
      (assert (< i (max-corr-length)))
      (let ([xs2 (take np (drop (- (max-corr-length) i) xs))])
        (jackknife
          + - *
          (map (lambda (xs) (apply average xs))
               (block-list (block-size) (map * xs1 xs2))))))))

(define (calc-corr xs)
  (let* ([corr (mk-corr xs)]
         [corr0 (corr 0)])
    (map
      (lambda (i)
        (tree-op / (corr i) corr0))
      (iota (max-corr-length)))))

(define (mk-change xs)
  (let* ([n (length xs)]
         [np (- n (max-corr-length))]
         [xs1 (drop (max-corr-length) xs)])
    (lambda (i)
      (let ([xs2 (take np (drop (- (max-corr-length) i) xs))])
        (map (lambda (xs) (apply average xs))
             (block-list (block-size) (map sqr (map - xs1 xs2))))))))

(define (calc-change xs)
  (let ([corr (mk-change xs)])
    (map corr (iota (max-change-length)))))

(define (mk-jackkife-list xs)
  (define (go i xs)
    (cond
      [(null? xs) (list)]
      [(<= i 0) (cdr xs)]
      [else (cons (car xs) (go (dec i) (cdr xs)))]))
  (let ([n (length xs)])
    (cons xs (map (lambda (i) (go i xs)) (iota n)))))

(define (avg-err vs)
  (list (apply average vs) (apply average-sigma vs)))

(define (j-avg-err js)
  (list (car js) (apply jackknife-sigma js)))

; -----------------------------------------

(define (get-force-path)
  "results/gm_force_info")

(define (get-force-fns)
  (list-sort
    (on < car)
    (filter
      (lambda (p)
        (and
          (>= (car p) (min-traj))
          (or (eq? #f (max-traj-num))
              (< (car p) (+ (min-traj) (max-traj-num))))))
      (map
        (lambda (p)
          (cons (string->number (car p)) (cdr p)))
        (paths->ppairs
          "traj="
          (glob-expand-at (get-force-path) "traj=*.lat")
          ".lat")))))

(define (get-force-data)
  (map
    (lambda (p)
      (load-lat-data (cdr p)))
    (get-force-fns)))

(define (get-force force-data)
  (let* ([_ (assert (not (null? force-data)))]
         [dim-sizes (lat-data-dim-sizes (car force-data))]
         [n-elem (vector-ref dim-sizes 1)])
    (map
      (lambda (i-elem)
        ((lambda (xs)
           (list i-elem
                 (apply average xs)
                 (apply std-deviation xs)))
         (apply append
           (map
             (lambda (ld)
               (let ([n-meas (vector-ref (lat-data-dim-sizes ld) 0)])
                 (map
                   (lambda (i-meas)
                     (vector-head (lat-data-ref ld (vector i-meas i-elem))))
                   (iota n-meas))))
             force-data))))
      (iota n-elem))))

; -----------------------------------------

(define (analysis-force)
  (let* ([force-data (get-force-data)]
         [force-list (get-force force-data)]
         )
    (print
      (list
        (list 'n (length force-data))
        (list 'force-list force-list)
        ))
    (mkdir-p "stats/info")
    (save-obj
      "stats/info/force.txt"
      (list
        (list 'n (length force-data))
        (list 'force-list force-list)
        ))))

(define (analysis)
  (let* ([data (get-data)]
         [b-data (block-list (block-size) (drop (max-corr-length) data))]
         [j-data (map (lambda (b) (apply append b)) (mk-jackkife-list b-data))]
         [avg (avg-err (map (lambda (d) (apply average d)) b-data))]
         [std-deviation
           (j-avg-err
             (map (lambda (vs) (sqrt (apply average (map sqr vs))))
                  (map subtract-avg j-data)))]
         [xs (subtract-avg data)]
         [corr (calc-corr xs)]
         [change (calc-change data)]
         )
    (print
      (list
        (list 'data-tag (data-tag))
        (list 'n (length data))
        (list 'n-block (length b-data))
        (list 'avg avg)
        (list 'std-deviation std-deviation)
        (list 'change (take 11 (map avg-err change)))
        (list 'partial-sum-corr (take 20 (map j-avg-err (partial-sum tree+ corr))))
        ))
    (mkdir-p "stats/info")
    (save-obj
      (format "stats/info/~a.txt" (data-tag))
      (list
        (list 'data-tag (data-tag))
        (list 'min-traj (min-traj))
        (list 'max-corr-length (max-corr-length))
        (list 'max-change-length (max-change-length))
        (list 'block-size (block-size))
        (list 'n (length data))
        (list 'n-block (length b-data))
        (list 'avg avg)
        (list 'std-deviation std-deviation)
        (list 'change (map avg-err change))
        (list 'corr (map j-avg-err corr))
        (list 'partial-sum-corr (map j-avg-err (partial-sum tree+ corr)))
        ; (list 'xs xs)
        (list 'all-data (get-all-data))
        ))))

(define (set-data-selector-plaq n-step)
  (data-tag (format "plaq-action-ape-smear-~a" n-step))
  (data-avg #f)
  (data-selector
    (lambda (ld)
      (vector-head
        (lat-data-ref ld (vector n-step 0))))))

(define (set-data-selector-topo n-step)
  (data-tag (format "topological-charge-ape-smear-~a" n-step))
  (data-avg 0)
  (data-selector
    (lambda (ld)
      (vector-head
        (lat-data-ref ld (vector n-step 2))))))

(define (all-analysis)
  (analysis-force)
  (set-data-selector-plaq 0)
  (analysis)
  (set-data-selector-plaq 1)
  (analysis)
  (set-data-selector-plaq 2)
  (analysis)
  (set-data-selector-topo 0)
  (analysis)
  (set-data-selector-topo 1)
  (analysis)
  (set-data-selector-topo 2)
  (analysis)
  )

(min-traj 2)

(max-traj-num #f)

(block-size 1)

(max-corr-length 1)

(max-change-length 1)

(all-analysis)
